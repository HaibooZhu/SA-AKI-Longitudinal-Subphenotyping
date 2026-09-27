# Original-label AutoGluon parameter comparison (2026-09-21)

User authorized remote training, Conda environment inspection, latest stable AutoGluon and parameter tuning. This experiment does not modify manuscript files or select new clustering labels.

Private server run directory: `~/jtim_autogluon_20260921` on the user-authorized SSH host. Original source directory: `~/fjmu_hanli_bak/kidney_sepsis_penotype_v3` (read-only throughout).

- `prepare.py` verifies original label hashes, labels against archived classifier inputs, cohort membership and disjoint train/test partitions. It constructs 243 predictors from 27 variables in time indices 1–4 (0–24 h), with max/min/mean and six pairwise differences. Results remain on the private server.
- Feature construction was independently compared against the supplied dynamic-window experiment's `shared_feature_builder.py` on all three cohorts. Largest numerical difference: 5.68e-14; labels identical after pointing the reference builder to original labels.
- `train.py` specifies 15 parameter configurations across XGBoost, LightGBM, CatBoost and ExtraTrees. All use the same original training population, five-fold bagging, one repeat, no stacking, seed 20260921, 16 CPU threads and no GPU. ExtraTrees' native out-of-bag shortcut is explicitly disabled. A weighted ensemble is also fitted from training out-of-fold probabilities.
- No supervised feature filtering is performed before cross-validation. All 243 raw predictors are supplied, unlike the older correlation/MI-filtered pipeline. Thus this is a new classifier experiment on the same phenotypes and observation window, not an exact one-parameter reproduction of the old classifier.
- Best XGBoost and best overall are locked by training out-of-fold macro OVO AUC before holdout evaluation. Original internal n=1,227 and external AUMC n=2,183 are never passed as tuning data. Training n=4,903.
- `verify_run.py` checks actual fold counts, independently recomputes selected out-of-fold and held-out metrics, and checks prediction labels and selection timing. Patient-level predictions and model files stay private; local reporting uses aggregate files only.

Environment: isolated Conda prefix under the run directory, Python 3.12, AutoGluon tabular 1.6.3 (latest stable PyPI release checked 2026-09-21; released 2026-09-18). Old base/main/jupy environments are preserved. Exact dependencies are saved in `environment.lock.txt`; `pip check` passed.

An initial training attempt was stopped before any held-out evaluation to disable ExtraTrees' default OOB shortcut and enforce the same five-fold design. Its script/log/partial outputs are preserved in `attempt_00/` on the server. No result-driven parameter change was made.

Official references checked during preparation:
- https://pypi.org/project/autogluon.tabular/1.6.3/
- https://auto.gluon.ai/stable/api/autogluon.tabular.TabularPredictor.fit.html

Current run results and status belong in private `02_revision_outputs/reports/20260921_autogluon_original_labels/`; consult the result report rather than treating this plan as proof of completion or improvement. No manuscripts, responses, public code or frozen inputs are changed by this experiment.
