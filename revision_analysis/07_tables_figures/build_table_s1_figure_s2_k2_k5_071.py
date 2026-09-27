#!/usr/bin/env python3
"""071: Table S1 and Figure S2 restricted to the candidate models present in the archived code (K = 2-5).

Source: the frozen mixAK notebooks (publication_code/0{1,2,3}.*/step2_mixAK_3-Fb2.py) fit mod2-mod5 only and
record their diagnostics as literal vectors (clusternum, mcmc_mus_failed, uncertainly_clustered_individuals,
deviance). This script parses those vectors, cross-checks eICU against the dynamic recomputation from the
archived model objects (reports/W3_archived_k_selection), and recomputes the composite Euclidean distance with
the archived rule: min-max rescaling of deviance and of the high-autocorrelation fraction, then the Euclidean norm.
No patient data are read. Outputs: Table_S1_k2_k5.csv and Figure_S2_k2_k5.{png,pdf,svg}.
"""
import argparse, math, re, json, hashlib
from pathlib import Path
import pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

EXP = Path(__file__).resolve().parents[2]
CODE = EXP / '00_frozen_inputs/code_snapshot/SA-AKI_Longitudinal_Subphenotype/publication_code'
SRC = {'MIMIC-IV': '01.MIMICIV_SAKI_trajCluster/step2_mixAK_3-Fb2.py',
       'AUMC': '02.AUMCdb_SAKI_trajCluster/step2_mixAK_3-Fb2.py',
       'eICU': '03.eICU_SAKI_trajCluster/step2_mixAK_3-Fb2.py'}
ARCH = EXP / '02_revision_outputs/reports/W3_archived_k_selection/archived_eicu_k2_k5_diagnostics.csv'


def vec(text, name):
    m = re.findall(r'^%s <- c\(([^)]*)\)' % re.escape(name), text, re.M)
    assert len(m) == 1, (name, len(m))
    return [float(v) for v in m[0].split(',')]


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    rows = []; prov = {}
    for co, rel in SRC.items():
        p = CODE / rel; t = p.read_text(encoding='utf-8')
        assert len(re.findall(r'mod(\d) <- GLMM_MCMC', t)) == 4 and 'mod6' not in t
        k = [int(v) for v in vec(t, 'clusternum')]; a = vec(t, 'mcmc_mus_failed')
        u = [int(v) for v in vec(t, 'uncertainly_clustered_individuals')]; d = vec(t, 'deviance')
        assert sorted(k) == [2, 3, 4, 5]
        dmin, dmax, amin, amax = min(d), max(d), min(a), max(a)
        for kk, aa, uu, dd in sorted(zip(k, a, u, d)):
            ed = math.hypot((dd - dmin) / (dmax - dmin), (aa - amin) / (amax - amin))
            rows.append({'cohort': co, 'K': kk, 'high_autocorr_pct': 100 * aa, 'uncertain_n': uu,
                         'deviance': dd, 'composite_ed': ed})
        prov[co] = {'file': str(p.relative_to(EXP)), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
    tab = pd.DataFrame(rows)
    # eICU cross-check with the archived model objects (dynamic recomputation, W3)
    arch = pd.read_csv(ARCH)
    kcol = next(c for c in arch.columns if c.lower() in ('k', 'clusternum', 'clusters'))
    e = tab[tab.cohort == 'eICU'].set_index('K')
    for r in arch.itertuples():
        kk = int(getattr(r, kcol))
        assert abs(e.loc[kk, 'deviance'] - round(float(r.mean_deviance), 1)) < 0.06, (kk, r.mean_deviance)
        assert int(e.loc[kk, 'uncertain_n']) == int(r.uncertain_patients)
        assert abs(e.loc[kk, 'high_autocorr_pct'] / 100 - float(r.high_absolute_lag1_fraction)) < 1e-9
    tab.to_csv(out / 'Table_S1_k2_k5.csv', index=False)
    (out / 'provenance.json').write_text(json.dumps({'sources': prov, 'eicu_crosscheck': str(ARCH.relative_to(EXP))}, indent=2))

    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 9,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.67))
    for ax, co in zip(axs, SRC):
        z = tab[tab.cohort == co]
        ax.plot(z.K, z.composite_ed, color='black', lw=1.6, marker='o', ms=6, mfc='white', mec='black', mew=1.6)
        ax.axvline(3, color='#E7211A', ls=(0, (4.5, 3)), lw=0.9)
        ax.text(3.08, 0.03, 'Optimal number of clusters: 3', fontsize=8.5, va='bottom')
        ax.set_title(co, fontsize=10.5)
        ax.set_xticks([2, 3, 4, 5]); ax.set_xlim(1.8, 5.2)
        ax.set_ylim(-0.05, 1.5); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
        ax.set_yticklabels([f'{v:.1f}' for v in [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]], rotation=90, va='center')
        ax.set_xlabel('Number of clusters'); ax.set_ylabel('Composite Euclidean distance')
        ax.tick_params(direction='out', length=3)
    fig.tight_layout(w_pad=2.0)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(out / f'Figure_S2_k2_k5.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(tab.round(3).to_string(index=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--out', type=Path, required=True)
    main(ap.parse_args().out)
