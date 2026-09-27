"""Reproduce frozen predictions and explain fold margins without densifying CSR.

No training or input mutation. Outputs go to a new, private audit directory.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '4')
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import xgboost as xgb
from autogluon.tabular import TabularPredictor


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    w = args.root / 'window_24h'
    out = w / 'shap_sparse_verified'
    out.mkdir(exist_ok=False)
    source_files = sorted((w / 'predictor').rglob('xgb.ubj')) + list(w.glob('private_predictions_*.pkl'))
    before = {str(p.relative_to(w)): digest(p) for p in source_files}
    predictor = TabularPredictor.load(str(w / 'predictor'))
    model = 'XGBoost_BAG_L1'
    bag = predictor._trainer.load_model(model)
    assert len(bag.models) == 5
    rep = {'xgboost': xgb.__version__, 'model': model, 'fold_count': 5, 'splits': {}, 'folds': []}
    contributions = []
    for split in ['internal', 'external']:
        data = pd.read_pickle(w / 'data' / f'{split}.pkl')
        saved = pd.read_pickle(w / f'private_predictions_{split}.pkl')
        assert np.array_equal(data.stay_id, saved.stay_id)
        X = data.drop(columns=['stay_id', 'dataset', 'groupHPD'])
        expected = saved[['p_DR', 'p_RR', 'p_PW']].to_numpy()
        predicted = predictor.predict_proba(X, model=model).loc[:, [1, 2, 3]].to_numpy()
        assert np.max(np.abs(predicted - expected)) < 1e-6
        Xt = predictor.transform_features(X)
        sparse_probs, dense_probs = [], []
        for name in bag.models:
            child = bag.load_child(name)
            child_prob = child.predict_proba(Xt.copy(), num_cpus=4)
            Xp = child.preprocess(Xt.copy())
            assert sp.isspmatrix_csr(Xp)
            features = list(child._features_internal or child.features)
            # No categorical expansion in this model: verify exact feature alignment.
            assert Xp.shape[1] == len(features)
            np.testing.assert_allclose(Xp.toarray(), Xt[features].to_numpy(), equal_nan=True)
            bst = child.model.get_booster()
            bst.set_param({'nthread': 4})
            iterations = (0, int(bst.best_iteration) + 1)
            dm = xgb.DMatrix(Xp)
            sparse_prob = bst.predict(dm, iteration_range=iterations)
            np.testing.assert_allclose(sparse_prob, child_prob, rtol=0, atol=1e-7)
            # Separately load the archived UBJ to establish file identity.
            file_bst = xgb.Booster()
            file_bst.load_model(str(Path(child.path) / 'xgb.ubj'))
            file_bst.set_param({'nthread': 4})
            np.testing.assert_allclose(file_bst.predict(dm, iteration_range=iterations), sparse_prob, rtol=0, atol=1e-7)
            dense_prob = bst.predict(xgb.DMatrix(Xp.toarray()), iteration_range=iterations)
            sparse_probs.append(sparse_prob)
            dense_probs.append(dense_prob)
            if split == 'internal':
                C = bst.predict(dm, pred_contribs=True, iteration_range=iterations)
                margin = bst.predict(dm, output_margin=True, iteration_range=iterations)
                error = float(np.max(np.abs(C.sum(axis=-1) - margin)))
                assert error < 5e-5, error
                contributions.append(C)
                rep['folds'].append({'name': name, 'best_iteration': iterations[1] - 1,
                                     'input_format': 'CSR', 'raw_booster_vs_child_max_abs_diff': float(np.max(np.abs(sparse_prob-child_prob))),
                                     'shap_margin_additivity_max_abs_diff': error})
        sparse_avg = np.mean(sparse_probs, axis=0)
        dense_avg = np.mean(dense_probs, axis=0)
        assert np.max(np.abs(sparse_avg - expected)) < 1e-6
        rep['splits'][split] = {'n': len(data), 'predictor_vs_saved_max_abs_diff': float(np.max(np.abs(predicted-expected))),
            'csr_fold_average_vs_saved_max_abs_diff': float(np.max(np.abs(sparse_avg-expected))),
            'dense_fold_average_vs_saved_max_abs_diff': float(np.max(np.abs(dense_avg-expected))),
            'dense_fold_average_vs_saved_mean_abs_diff': float(np.mean(np.abs(dense_avg-expected)))}
        print(split, rep['splits'][split], flush=True)
    C = np.stack(contributions)
    np.save(out / 'shap_contribs_internal_5fold_PRIVATE.npy', C)
    mean_abs = np.abs(C[..., :-1]).mean(axis=1)
    importance = pd.DataFrame({'feature': features, **{label: mean_abs[:, i].mean(0) for i, label in enumerate(['DR', 'RR', 'PW'])}})
    importance['total'] = importance[['DR', 'RR', 'PW']].sum(axis=1)
    importance['fold_sd_total'] = mean_abs.sum(axis=1).std(axis=0)
    importance = importance.sort_values('total', ascending=False).reset_index(drop=True)
    importance['rank'] = importance.index + 1
    importance.to_csv(out / 'shap_mean_abs_by_feature_24h.csv', index=False)
    rep['shap_scale'] = 'Per-class raw margin; mean absolute SHAP across internal-test patients, then averaged over five folds. Not a decomposition of bagged probabilities.'
    rep['cause'] = 'Direct dense input bypasses the AutoGluon child CSR transformation; implicit CSR zeros and explicitly stored dense zeros have different missing-value semantics.'
    rep['source_hashes'] = before
    rep['sources_unchanged'] = all(digest(w / name) == sha for name, sha in before.items())
    assert rep['sources_unchanged']
    (out / 'shap_verification.json').write_text(json.dumps(rep, indent=2))
    (out / 'COMPLETE.json').write_text(json.dumps({'verified': True, 'n_features': len(features), 'contribution_shape': list(C.shape)}))
    print(importance.head(10).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
