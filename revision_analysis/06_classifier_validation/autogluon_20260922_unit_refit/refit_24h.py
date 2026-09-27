"""Unit-harmonized 24 h refit; fixed original labels, splits and XGB settings.

Private inputs/predictions/models stay inside the authorized server run folder.
This is a paired processing comparison, not external-test model selection.
"""
import argparse
import hashlib
import importlib.metadata as md
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from autogluon.tabular import TabularPredictor
from sklearn.metrics import roc_auc_score


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str) + '\n')


def normalize(tr, cohort):
    out = tr.copy()
    h = pd.to_numeric(out.hematocrit, errors='coerce')
    fraction = h.between(.05, .80)
    pct = h.between(5, 80)
    out['hematocrit'] = h.where(pct)
    out.loc[fraction, 'hematocrit'] = h.loc[fraction] * 100
    f = pd.to_numeric(out.fio2, errors='coerce')
    if cohort == 'eicu':
        assert not f.dropna().gt(1).any(), 'Unexpected eICU FiO2 scale'
        f = f * 100
    out['fio2'] = f.where(f.between(21, 100))
    # Do not fill AUMC bilirubin from the post-imputation trajectory matrix.
    return out


def summary(series):
    x = pd.to_numeric(series, errors='coerce')
    return {'nonmissing': int(x.notna().sum()), 'missing': int(x.isna().sum()),
            'median': None if x.dropna().empty else float(x.median()),
            'min': None if x.dropna().empty else float(x.min()),
            'max': None if x.dropna().empty else float(x.max())}


def main(args):
    os.umask(0o077)
    sys.path.insert(0, str(args.legacy_scripts))
    import run_windows as rw
    from prepare import sha
    from train import metrics

    root = args.root
    root.mkdir(parents=True, exist_ok=False)
    (root / 'data').mkdir()
    seed = 20260921
    plan = {'window_hours': 24, 'parameters': rw.HP, 'seed': seed,
            'num_bag_folds': 5, 'num_bag_sets': 1, 'num_stack_levels': 0,
            'num_cpus': 16, 'training_n': 4903, 'internal_n': 1227, 'external_n': 2183,
            'changes': ['eICU FiO2 fraction to percent before feature aggregation',
                        'all-cohort FiO2 values outside 21-100 set to missing',
                        'hematocrit 0.05-0.80 times 100; 5-80 retained; others missing'],
            'bilirubin': 'No new imputation or field replacement; absent in AUMC pre-MI archive',
            'selection': 'One fixed-parameter refit; no selection using test performance',
            'unit_rules_source': 'Existing regenerate_supplementary_tables.py',
            'baseline_run': str(args.baseline)}
    dump(root / 'plan.json', plan)
    dump(root / 'environment.json', {'python': platform.python_version(),
        'versions': {p: md.version(p) for p in ['autogluon.tabular', 'autogluon.core',
                     'xgboost-cpu', 'numpy', 'pandas', 'scikit-learn']}})
    manifest = json.loads((args.baseline / 'input_manifest.json').read_text())
    for path, expected in manifest['sources'].items():
        assert sha(Path(path)) == expected, f'Source changed: {path}'
    audit = {'original_source_hashes_verified': True, 'inputs': manifest['sources'],
             'processing': {}, 'splits': {}, 'temporal_boundary_checks': []}
    original_frames, corrected_frames = [], []
    for co, (folder, oldtraj, labfile, n) in rw.SPECS.items():
        lab = pd.read_csv(args.source / folder / labfile, usecols=['stay_id', 'groupHPD']).astype('int64')
        assert lab.groupby('stay_id').groupHPD.nunique().max() == 1
        lab = lab.drop_duplicates('stay_id').set_index('stay_id')
        assert len(lab) == n
        tr = pd.read_csv(args.source / folder / 'df_im_By_ffill.csv')
        tr['stay_id'] = tr.stay_id.astype('int64')
        if co == 'mimic':
            tr = tr.rename(columns={c: c[:-5] for c in tr if c.endswith('_mean')})
        if 'bilirubin' not in tr and 'bilirubin_total' in tr:
            tr['bilirubin'] = tr.bilirubin_total
        bili_present = 'bilirubin' in tr
        bp, idcol, bcol, scale = rw.BASELINE[co]
        baseline = pd.read_csv(args.source / bp, usecols=[idcol, bcol]).rename(columns={idcol: 'stay_id', bcol: 'baseline'})
        baseline['stay_id'] = baseline.stay_id.astype('int64')
        baseline['baseline'] *= scale
        tr = tr.merge(baseline, on='stay_id', how='left', validate='many_to_one')
        assert tr.loc[tr.stay_id.isin(lab.index), 'baseline'].notna().all()
        tr['crea_divide_basecrea'] = (tr.creatinine / tr.baseline).round(2)
        for f in rw.BASE:
            tr[f] = pd.to_numeric(tr[f], errors='coerce') if f in tr else np.nan
        tr = tr[['stay_id', 'time'] + rw.BASE]
        fixed = normalize(tr, co)
        allowed = ['hematocrit', 'fio2']
        pd.testing.assert_frame_equal(tr.drop(columns=allowed), fixed.drop(columns=allowed))
        scope = tr.stay_id.isin(lab.index) & tr.time.between(1, 4)
        audit['processing'][co] = {'bilirubin_column_present': bili_present, 'fields': {}}
        for f in allowed + ['bilirubin']:
            audit['processing'][co]['fields'][f] = {
                'before': summary(tr.loc[scope, f]), 'after': summary(fixed.loc[scope, f]),
                'changed_rows': int((~np.isclose(tr.loc[scope, f], fixed.loc[scope, f], equal_nan=True)).sum())}
        original = rw.features(tr, lab.index, 4)
        corrected = rw.features(fixed, lab.index, 4)
        perturbed = tr.copy()
        perturbed.loc[perturbed.time > 4, rw.BASE] = 987654.321
        # Normalize only the eligible window for the future-perturbation test.
        # Values beyond 24 h are never consumed, including by the scale assertion.
        perturbed = perturbed.loc[perturbed.time <= 4].copy()
        pd.testing.assert_frame_equal(corrected, rw.features(normalize(perturbed, co), lab.index, 4))
        audit['temporal_boundary_checks'].append({'cohort': co, 'passed': True})
        for frame, collection in [(original, original_frames), (corrected, corrected_frames)]:
            frame['groupHPD'] = lab.groupHPD
            frame['dataset'] = co
            collection.append(frame)
    original_all = pd.concat(original_frames)
    corrected_all = pd.concat(corrected_frames)
    assert original_all.index.is_unique and corrected_all.index.is_unique
    for split, n in [('train', 4903), ('internal', 1227), ('external', 2183)]:
        prior = pd.read_pickle(args.baseline / 'window_24h/data' / f'{split}.pkl')
        old = original_all.loc[prior.stay_id].reset_index()
        pd.testing.assert_frame_equal(old[prior.columns], prior, check_dtype=False)
        fixed = corrected_all.loc[prior.stay_id].reset_index()[prior.columns]
        assert len(fixed) == n
        pd.testing.assert_frame_equal(fixed[['stay_id', 'dataset', 'groupHPD']], prior[['stay_id', 'dataset', 'groupHPD']])
        changed = []
        for c in prior.columns:
            if c in ['stay_id', 'dataset', 'groupHPD']:
                continue
            if not np.allclose(prior[c], fixed[c], equal_nan=True):
                changed.append(c)
                assert c.startswith(('hematocrit_', 'fio2_')), c
        fixed.to_pickle(root / 'data' / f'{split}.pkl')
        audit['splits'][split] = {'n': n, 'original_features_reproduced': True,
            'members_and_labels_unchanged': True, 'changed_columns': changed,
            'class_counts': fixed.groupHPD.value_counts().sort_index().to_dict(),
            'file_sha256': sha(root / 'data' / f'{split}.pkl')}
    dump(root / 'input_verification.json', audit)
    print('PREPARATION PASSED: original features reproduced; only Hct/FiO2 changed', flush=True)
    tr = pd.read_pickle(root / 'data/train.pkl')
    x = tr.drop(columns=['stay_id', 'dataset'])
    predictor = TabularPredictor(label='groupHPD', problem_type='multiclass',
        eval_metric='roc_auc_ovo_macro', path=str(root / 'predictor'), verbosity=2,
        learner_kwargs={'random_state': seed})
    predictor.fit(x, hyperparameters={'XGB': rw.HP}, num_bag_folds=5, num_bag_sets=1,
        num_stack_levels=0, num_cpus=16, num_gpus=0,
        ag_args_ensemble={'fold_fitting_strategy': 'sequential_local'},
        fit_weighted_ensemble=False, calibrate=False, calibrate_decision_threshold=False,
        time_limit=1200)
    lb = predictor.leaderboard(silent=True)
    assert len(lb) == 1
    lb.to_csv(root / 'training_cv.csv', index=False)
    model = lb.iloc[0].model
    bag = predictor._trainer.load_model(model)
    assert len(bag.models) == 5
    folds = list(bag._cv_splitters[0].split(x.drop(columns='groupHPD'), tr.groupHPD))
    digest = hashlib.sha256(b''.join(np.asarray(i, dtype='int64').tobytes() + np.asarray(j, dtype='int64').tobytes() for i, j in folds)).hexdigest()
    original_folds = json.loads((args.baseline / 'window_24h/fit_details.json').read_text())['fold_hash']
    assert digest == original_folds
    oof = predictor.predict_proba_oof(model=model).loc[:, [1, 2, 3]].to_numpy()
    assert abs(metrics(tr.groupHPD.to_numpy(), oof)['auc_ovo_macro'] - float(lb.iloc[0].score_val)) < 1e-10
    dump(root / 'fit_verification.json', {'model': model, 'five_fold_hash': digest,
        'matches_baseline_folds': True, 'oof_auc_recomputed': True,
        'hyperparameters': predictor.model_hyperparameters(model)})
    rows, predictions = [], {}
    for split in ['internal', 'external']:
        test = pd.read_pickle(root / 'data' / f'{split}.pkl')
        y = test.groupHPD.to_numpy(dtype=int)
        features = test.drop(columns=['stay_id', 'dataset', 'groupHPD'])
        probs = predictor.predict_proba(features, model=model).loc[:, [1, 2, 3]].to_numpy()
        assert np.isfinite(probs).all() and np.allclose(probs.sum(axis=1), 1, atol=1e-5)
        saved = pd.DataFrame({'stay_id': test.stay_id, 'groupHPD': y,
                            'p_DR': probs[:, 0], 'p_RR': probs[:, 1], 'p_PW': probs[:, 2]})
        saved.to_pickle(root / f'private_predictions_{split}.pkl')
        old = pd.read_pickle(args.baseline / 'window_24h' / f'private_predictions_{split}.pkl')
        pd.testing.assert_frame_equal(old[['stay_id', 'groupHPD']], saved[['stay_id', 'groupHPD']])
        oldp = old[['p_DR', 'p_RR', 'p_PW']].to_numpy()
        for version, p in [('original', oldp), ('unit_harmonized', probs)]:
            rows.append({'split': split, 'version': version, 'n': len(y), **metrics(y, p)})
        predictions[split] = (y, oldp, probs)
    pd.DataFrame(rows).to_csv(root / 'comparison_metrics.csv', index=False)
    # Verify probabilities using the saved AutoGluon predictor (including its CSR preprocessing).
    loaded = TabularPredictor.load(str(root / 'predictor'))
    for split in predictions:
        test = pd.read_pickle(root / 'data' / f'{split}.pkl')
        actual = loaded.predict_proba(test.drop(columns=['stay_id', 'dataset', 'groupHPD']), model=model).loc[:, [1, 2, 3]].to_numpy()
        assert np.max(np.abs(actual - predictions[split][2])) < 1e-7
    dump(root / 'COMPLETE.json', {'evaluations': len(rows), 'one_fixed_refit': True,
        'saved_predictor_reproduces_probabilities': True, 'script_sha256': sha(Path(__file__))})
    print(pd.DataFrame(rows)[['split', 'version', 'auc_ovo_macro', 'auc_DR_RR', 'auc_DR_PW', 'auc_RR_PW']].to_string(index=False), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['source', 'baseline', 'legacy-scripts', 'root']:
        p.add_argument('--' + name, type=Path, required=True)
    main(p.parse_args())
