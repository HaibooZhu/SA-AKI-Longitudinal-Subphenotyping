"""Compare K = 2 and K = 3 solutions from the recovered W21 extended-sampling refits (no refitting).

Supports the note to Supplementary Table S17A. For each cohort it reports, as aggregate counts only:
1. how the two K = 2 clusters map onto the original subphenotypes (RR/DR/PW; originally uncertain patients separately);
2. agreement of each K = 3 initialization with the original assignments;
3. the posterior mean deviance of K = 2 minus that of the K = 3 initializations reproducing the original assignments,
   compared with a Bayesian-information-criterion-type penalty for the additional component parameters
   (number of additional parameters x ln[number of observations]). The deviance is the Laplace-approximated
   marginal deviance reported by mixAK::GLMM_MCMC.
Patient-level inputs are read in the authorized environment and never written.
"""
from pathlib import Path
import argparse, itertools, json
import numpy as np, pandas as pd

MATRICES = {'mimic': '01.MIMICIV_SAKI_trajCluster/df_mixAK_fea4_C3.csv',
            'aumc': '02.AUMCdb_SAKI_trajCluster/df_mixAK_fea3_C3_aumc.csv',
            'eicu': '03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv'}
Q_DIM = {'mimic': 8, 'aumc': 6, 'eicu': 8}          # random-effect dimension = indicators x (intercept + slope)
SEEDS = [20260805, 20260806, 20260807]


def align(ref, cand, k):
    cross = pd.crosstab(cand, ref).reindex(index=range(1, k + 1), columns=range(1, k + 1), fill_value=0)
    best = max(itertools.permutations(range(k)), key=lambda p: sum(cross.iloc[i, p[i]] for i in range(k)))
    return {i + 1: best[i] + 1 for i in range(k)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work-root', type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument('--output', type=Path, default=None)
    a = ap.parse_args()
    rec = a.work_root / '02_revision_outputs/reports/W21_burn_in_remediation/recovered_final_20260908/outputs'
    snap = a.work_root / '00_frozen_inputs/data_snapshot/remote_project_snapshot'
    output = a.output or a.work_root / '02_revision_outputs/reports/W21_burn_in_remediation/k2_k3_crosstab_summary.json'
    out = {}
    for co, path in MATRICES.items():
        matrix = pd.read_csv(snap / path, usecols=['stay_id', 'groupHPD'])
        orig = matrix.drop_duplicates('stay_id').copy()
        prev = orig[orig.groupHPD.isin([1, 2, 3])].groupHPD.value_counts().sort_values(ascending=False)
        name = {prev.index[0]: 'RR', prev.index[1]: 'DR', prev.index[2]: 'PW'}
        orig['orig'] = orig.groupHPD.map(name).fillna('uncertain')
        res = {'n_patients': len(orig), 'n_observations': len(matrix), 'orig_counts': orig.orig.value_counts().to_dict()}

        dev = {(k, s): float(pd.read_csv(rec / f'{co}__primary_converged__K{k}__seed{s}__diagnostics.csv').iloc[0].mean_deviance)
               for k in (2, 3) for s in SEEDS}

        # 1. K = 2 clusters versus the original subphenotypes
        k2 = {}
        for s in SEEDS:
            x = pd.read_csv(rec / f'{co}__primary_converged__K2__seed{s}__assignments.csv').merge(orig, on='stay_id', validate='one_to_one')
            assert len(x) == len(orig)
            big = x.group_median.value_counts().idxmax()
            x['K2'] = np.where(x.group_median == big, 'K2_large', 'K2_small')
            k2[s] = x.set_index('stay_id').K2
            res[f'K2_seed{s}_x_orig'] = pd.crosstab(x.K2, x.orig).to_dict()
        res['K2_seed_pairwise_agreement'] = [float((k2[SEEDS[0]] == k2[s]).mean()) for s in SEEDS[1:]]

        # 2. each K = 3 initialization versus the original assignments
        for s in SEEDS:
            x = pd.read_csv(rec / f'{co}__primary_converged__K3__seed{s}__assignments.csv').merge(orig, on='stay_id', validate='one_to_one')
            mp = align(x.groupHPD, x.group_median, 3)
            res[f'K3_seed{s}'] = {'agreement_vs_original': float((x.group_median.map(mp) == x.groupHPD).mean()),
                                  'mean_deviance': dev[(3, s)]}

        # 3. deviance difference versus the penalty for one additional component
        q = Q_DIM[co]; dp = q + q * (q + 1) // 2 + 1        # mean q, covariance q(q+1)/2, weight 1
        reproducing = [dev[(3, s)] for s in SEEDS if res[f'K3_seed{s}']['agreement_vs_original'] > 0.99]
        d2 = min(dev[(2, s)] for s in SEEDS)
        res['penalty_check'] = {'extra_params_per_component': dp, 'min_mean_dev_K2': d2,
                                'K3_reproducing_initializations': len(reproducing),
                                'min_dev_K2_minus_max_dev_reproducing_K3': d2 - max(reproducing),
                                'BIC_penalty_ln_observations': round(float(dp * np.log(len(matrix))), 1)}
        out[co] = res
    output.write_text(json.dumps(out, indent=1, default=str), encoding='utf-8')
    for co, r in out.items():
        pw = {s: r[f'K2_seed{s}_x_orig'].get('PW', {}).get('K2_small', 0) for s in SEEDS}
        print(co, 'patients', r['n_patients'], 'observations', r['n_observations'], 'PW in small K2 cluster', pw, r['penalty_check'])


if __name__ == '__main__':
    main()
