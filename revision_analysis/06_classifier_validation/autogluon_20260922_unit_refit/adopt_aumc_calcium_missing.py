#!/usr/bin/env python3
"""Adopt the AUMC calcium-missing external predictions (author decision 2026-09-23).

Inputs (server, read only):
  - fitted six-window models and saved predictions: ~/jtim_refit_units_20260922/windows
  - external predictions with every calcium_ feature set to missing (audit_aumc_calcium.py):
    ~/jtim_calcium_audit_20260923
No model is fitted or changed. Internal predictions are linked unchanged.

Steps:
  1. Rebuild a six-window evaluation root whose external predictions are the calcium-missing ones.
  2. Rerun the original bootstrap (analyze_windows.analyze; seed 20260922, 1000 class-stratified draws).
  3. 24 h evaluation as in evaluate_24h.py: pair-normalized AUCs and ROC points, class-specific
     precision/recall/F1 and one-vs-rest calibration intercept/slope (unpenalized logistic regression
     of the outcome on logit(p), fitted by IRLS because statsmodels is not installed on the server),
     for both the saved and the calcium-missing external predictions and for the internal set.
Aggregate CSV/JSON only are written to OUT/aggregate; patient-level files stay in OUT/windows.
"""
import os, sys, json, shutil, hashlib, time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_fscore_support, accuracy_score

SRC = Path('~/jtim_refit_units_20260922/windows').expanduser()
AUD = Path('~/jtim_calcium_audit_20260923').expanduser()
LEG = Path('~/jtim_autogluon_20260921').expanduser()
OUT = Path('~/jtim_calcium_adopted_20260923').expanduser()
HOURS = [6, 12, 18, 24, 30, 36]
CLASSES = [1, 2, 3]; NAME = {1: 'DR', 2: 'RR', 3: 'PW'}
PAIRS = [(1, 2, 'DR vs RR'), (1, 3, 'DR vs PW'), (2, 3, 'RR vs PW')]


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''): h.update(c)
    return h.hexdigest()


def calibration_parameters(yb, prob, iters=100):
    """Unpenalized logistic regression yb ~ a + b*logit(p) by IRLS (same model as statsmodels GLM Binomial)."""
    p = np.clip(prob, 1e-6, 1 - 1e-6)
    X = np.column_stack([np.ones_like(p), np.log(p / (1 - p))])
    beta = np.zeros(2)
    for _ in range(iters):
        eta = X @ beta; mu = 1 / (1 + np.exp(-eta)); w = mu * (1 - mu)
        step = np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (yb - mu))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-12: break
    return float(beta[0]), float(beta[1])


def evaluate(split, y, P, tag):
    out = {'summary': {'split': split, 'version': tag, 'n': int(len(y)),
                       'macro_ovo_auc': float(roc_auc_score(y, P, multi_class='ovo', average='macro', labels=CLASSES))},
           'pairs': [], 'roc': [], 'classes': []}
    for a, b, lab in PAIRS:
        m = np.isin(y, [a, b]); pa = P[m, a - 1]; pb = P[m, b - 1]; s = pa / np.maximum(pa + pb, 1e-15)
        out['pairs'].append({'split': split, 'version': tag, 'pair': lab, 'n': int(m.sum()),
                             'auc_pair_normalized': float(roc_auc_score(y[m] == a, s))})
        fpr, tpr, _ = roc_curve(y[m] == a, s)
        out['roc'] += [{'split': split, 'version': tag, 'pair': lab, 'fpr': float(f), 'tpr': float(t)} for f, t in zip(fpr, tpr)]
    pred = np.array(CLASSES)[P.argmax(1)]
    pr, rc, f1, sup = precision_recall_fscore_support(y, pred, labels=CLASSES, zero_division=0)
    out['summary']['accuracy'] = float(accuracy_score(y, pred))
    for i, c in enumerate(CLASSES):
        yb = (y == c).astype(float); ic, sl = calibration_parameters(yb, P[:, i])
        out['classes'].append({'split': split, 'version': tag, 'class': NAME[c], 'n': int(sup[i]),
                               'precision': float(pr[i]), 'recall': float(rc[i]), 'f1': float(f1[i]),
                               'auc_ovr': float(roc_auc_score(yb, P[:, i])),
                               'calibration_intercept': ic, 'calibration_slope': sl})
    return out


def main():
    os.umask(0o077)
    OUT.mkdir(exist_ok=False)
    root = OUT / 'windows'; root.mkdir()
    agg = OUT / 'aggregate'; agg.mkdir()
    sys.path.insert(0, str(LEG))
    from train import metrics
    import analyze_windows as aw

    watched = [SRC / 'metrics.csv', SRC / 'auc_confidence_intervals.csv', SRC / 'paired_changes_vs24h.csv']
    old = pd.read_csv(SRC / 'metrics.csv')
    rows = []; checks = {}
    for h in HOURS:
        s = (SRC / f'window_{h:02d}h').resolve()
        w = root / f'window_{h:02d}h'; w.mkdir()
        (w / 'data').symlink_to(s / 'data', target_is_directory=True)
        shutil.copy2(s / 'model_locked.json', w / 'model_locked.json')
        (w / 'private_predictions_internal.pkl').symlink_to(s / 'private_predictions_internal.pkl')
        watched += [s / 'private_predictions_internal.pkl', s / 'private_predictions_external.pkl',
                    AUD / f'private_predictions_external_{h:02d}h_calcium_missing.pkl']
        ext_old = pd.read_pickle(s / 'private_predictions_external.pkl')
        ext_new = pd.read_pickle(AUD / f'private_predictions_external_{h:02d}h_calcium_missing.pkl')
        pd.testing.assert_frame_equal(ext_old[['stay_id', 'groupHPD']].reset_index(drop=True),
                                      ext_new[['stay_id', 'groupHPD']].reset_index(drop=True))
        ext_new.to_pickle(w / 'private_predictions_external.pkl')
        for split in ['internal', 'external']:
            r = old[(old.split == split) & (old.hours == h)].iloc[0].to_dict()
            if split == 'external':
                y = ext_new.groupHPD.to_numpy(dtype=int); P = ext_new[['p_DR', 'p_RR', 'p_PW']].to_numpy()
                r.update(metrics(y, P))
            rows.append(r)
        checks[f'{h}h'] = {'external_ids_labels_order_unchanged': True,
                           'max_abs_prob_change_external': float(np.max(np.abs(
                               ext_old[['p_DR', 'p_RR', 'p_PW']].to_numpy() - ext_new[['p_DR', 'p_RR', 'p_PW']].to_numpy())))}
    before = {str(p): sha(p) for p in watched}
    pd.DataFrame(rows)[old.columns].to_csv(root / 'metrics.csv', index=False)
    (root / 'TRAINING_COMPLETE.json').write_text(json.dumps({'time': time.time(), 'windows': 6,
        'note': 'No training. Fitted models of jtim_refit_units_20260922 reused; external predictions from '
                'jtim_calcium_audit_20260923 (all calcium_ features set to missing in AUMC).'}))
    aw.analyze(root, 1000)

    # bootstrap outputs: internal must equal the previous run exactly (same seed, same predictions)
    ci_new = pd.read_csv(root / 'auc_confidence_intervals.csv'); ci_old = pd.read_csv(SRC / 'auc_confidence_intervals.csv')
    ch_new = pd.read_csv(root / 'paired_changes_vs24h.csv'); ch_old = pd.read_csv(SRC / 'paired_changes_vs24h.csv')
    for new_, old_, name in [(ci_new, ci_old, 'ci'), (ch_new, ch_old, 'paired')]:
        a = new_[new_.split == 'internal'].reset_index(drop=True); b = old_[old_.split == 'internal'].reset_index(drop=True)
        pd.testing.assert_frame_equal(a, b, check_exact=True)
        checks[f'internal_{name}_identical_to_previous_run'] = True
    m_new = pd.read_csv(root / 'metrics.csv')
    pd.testing.assert_frame_equal(m_new[m_new.split == 'internal'].reset_index(drop=True),
                                  old[old.split == 'internal'].reset_index(drop=True), check_exact=True)
    checks['internal_metrics_identical'] = True

    # 24 h evaluation
    s24 = (SRC / 'window_24h').resolve()
    res = []
    for split, path, tag in [('internal', s24 / 'private_predictions_internal.pkl', 'saved'),
                             ('external', s24 / 'private_predictions_external.pkl', 'saved'),
                             ('external', root / 'window_24h' / 'private_predictions_external.pkl', 'calcium_missing')]:
        d = pd.read_pickle(path)
        res.append(evaluate(split, d.groupHPD.to_numpy(dtype=int), d[['p_DR', 'p_RR', 'p_PW']].to_numpy(float), tag))
    pd.DataFrame([r['summary'] for r in res]).to_csv(agg / 'summary_24h.csv', index=False)
    pd.DataFrame(sum([r['pairs'] for r in res], [])).to_csv(agg / 'pairwise_auc_24h.csv', index=False)
    pd.DataFrame(sum([r['classes'] for r in res], [])).to_csv(agg / 'class_metrics_calibration_24h.csv', index=False)
    pd.DataFrame(sum([r['roc'] for r in res], [])).to_csv(agg / 'roc_curves_24h.csv', index=False)
    # cross-check 24 h evaluation against bootstrap point estimates
    for r in res:
        if r['summary']['version'] == 'saved' and r['summary']['split'] == 'external':
            continue
        z = ci_new[(ci_new.split == r['summary']['split']) & (ci_new.hours == 24)].set_index('metric').auc
        assert abs(z['auc_ovo_macro'] - r['summary']['macro_ovo_auc']) < 1e-9
        for pr_, key in zip(r['pairs'], ['auc_DR_RR', 'auc_DR_PW', 'auc_RR_PW']):
            assert abs(z[key] - pr_['auc_pair_normalized']) < 1e-9
    checks['24h_evaluation_matches_bootstrap_points'] = True

    for f in ['metrics.csv', 'auc_confidence_intervals.csv', 'paired_changes_vs24h.csv', 'evaluation_verification.json']:
        shutil.copy2(root / f, agg / f)
    after = {str(p): sha(p) for p in watched}
    assert before == after, 'source files changed'
    checks['source_files_unchanged'] = True
    checks['script_sha256'] = sha(Path(__file__))
    (agg / 'verification.json').write_text(json.dumps(checks, indent=2))
    (agg / 'COMPLETE.json').write_text(json.dumps({'no_fit_called': True, 'bootstrap_seed': 20260922,
                                                    'draws': 1000, 'time': time.time()}, indent=2))
    print(json.dumps(checks, indent=1))


if __name__ == '__main__':
    main()
