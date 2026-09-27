"""Complete five fixed-window refits, reuse verified 24 h, then evaluate all six."""
from pathlib import Path
import json
import subprocess
import sys
import time
import pandas as pd

ROOT = Path('~/jtim_refit_units_20260922').expanduser()
WINDOWS = ROOT / 'windows'
LEGACY = Path('~/jtim_autogluon_20260921').expanduser()
SOURCE = Path('~/fjmu_hanli_bak/kidney_sepsis_penotype_v3').expanduser()
BASELINE = Path('~/jtim_windows_20260921').expanduser()
WINDOWS.mkdir(exist_ok=False)
(WINDOWS / 'window_24h').symlink_to(ROOT / 'run_24h', target_is_directory=True)
for hours in [6, 12, 18, 30, 36]:
    subprocess.run([sys.executable, str(ROOT / 'code/refit_window.py'),
        '--source', str(SOURCE), '--baseline', str(BASELINE),
        '--legacy-scripts', str(LEGACY), '--root', str(WINDOWS / f'window_{hours:02d}h'),
        '--hours', str(hours)], check=True)

rows = []
folds = []
for hours in [6, 12, 18, 24, 30, 36]:
    w = WINDOWS / f'window_{hours:02d}h'
    assert json.loads((w / 'COMPLETE.json').read_text())['saved_predictor_reproduces_probabilities']
    fit = json.loads((w / 'fit_verification.json').read_text())
    folds.append(fit['five_fold_hash'])
    # One fixed model was specified in plan.json before training/evaluation.
    (w / 'model_locked.json').write_text(json.dumps({'model': fit['model'],
        'time': (w / 'plan.json').stat().st_mtime,
        'source': 'fixed-parameter plan written before fitting; no model selection'}))
    df = pd.read_csv(w / 'comparison_metrics.csv')
    df = df[df.version.eq('unit_harmonized')].drop(columns='version')
    df['hours'] = hours
    df['model'] = fit['model']
    rows.append(df)
assert len(set(folds)) == 1
pd.concat(rows, ignore_index=True).to_csv(WINDOWS / 'metrics.csv', index=False)
(WINDOWS / 'TRAINING_COMPLETE.json').write_text(json.dumps({'time': time.time(),
    'windows': 6, 'same_folds': True, '24h_reused_without_refitting': True}))
subprocess.run([sys.executable, str(LEGACY / 'analyze_windows.py'), '--root', str(WINDOWS)], check=True)
subprocess.run([sys.executable, str(ROOT / 'code/verify_sparse_shap_24h.py'), '--root', str(WINDOWS)], check=True)
