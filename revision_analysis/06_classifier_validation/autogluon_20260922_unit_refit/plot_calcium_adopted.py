#!/usr/bin/env python3
"""Figure 5, Figure S8 and Table S19/S20 CSVs from aggregate outputs of adopt_aumc_calcium_missing.py.

Plotting code is copied unchanged from evaluate_24h.py (Figure 5) and plot_windows.py (Figure S8);
only the data source differs. With --version saved and the 20260922 inputs it redraws the previous
figures (used to check that the drawing is unchanged).
Usage: plot_calcium_adopted.py --src <aggregate dir> --out <dir> [--external-version calcium_missing|saved]
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument('--src', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--external-version', default='calcium_missing')
p.add_argument('--old-layout', action='store_true', help='read 20260922 file layout (no version column)')
a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=True)

COHORT = {'internal': 'MIMIC-IV + eICU internal test set', 'external': 'AUMC external test set'}
PAIRS = [(1, 2, 'DR vs RR'), (1, 3, 'DR vs PW'), (2, 3, 'RR vs PW')]
if a.old_layout:
    ev = a.src / '24h_model_evaluation'
    rocs = pd.read_csv(ev / 'roc_curves_24h.csv'); pairs = pd.read_csv(ev / 'pairwise_auc_24h.csv')
    import json; summ = json.load(open(ev / 'evaluation_24h.json'))['summary']
    summary = {s: {'n': summ[s]['n'], 'macro_ovo_auc': summ[s]['macro_ovo_auc']} for s in summ}
else:
    keep = lambda d: d[(d.split == 'internal') | (d.version == a.external_version)]
    rocs = keep(pd.read_csv(a.src / 'roc_curves_24h.csv')); pairs = keep(pd.read_csv(a.src / 'pairwise_auc_24h.csv'))
    s = keep(pd.read_csv(a.src / 'summary_24h.csv'))
    summary = {r.split: {'n': int(r.n), 'macro_ovo_auc': float(r.macro_ovo_auc)} for r in s.itertuples()}

# ---------------- Figure 5 (code from evaluate_24h.py)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(11,5.2))
colors={'DR vs RR':'#1f77b4','DR vs PW':'#d62728','RR vs PW':'#2ca02c'}
for ax,split,letter in zip(axes,['internal','external'],['A','B']):
    for _a,_b,lab in PAIRS:
        sub=rocs[(rocs.split==split)&(rocs.pair==lab)]; auc=pairs[(pairs.split==split)&(pairs.pair==lab)].auc_pair_normalized.iloc[0]
        ax.plot(sub.fpr,sub.tpr,color=colors[lab],lw=2,label=f'{lab} (AUC {auc:.3f})')
    ax.plot([0,1],[0,1],'k--',lw=1.2)
    macro=summary[split]['macro_ovo_auc']
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate'); ax.set_xlim(-0.01,1.01); ax.set_ylim(-0.01,1.01)
    ax.set_title(f'{COHORT[split]} (n = {summary[split]["n"]:,})\nmacro one-versus-one AUC {macro:.3f}',fontsize=11)
    ax.legend(loc='lower right',frameon=False,title='Subphenotype pair'); ax.set_aspect('equal')
    ax.text(-0.12,1.06,letter,transform=ax.transAxes,fontsize=18,fontweight='bold')
fig.tight_layout()
for ext in ['pdf','svg','png']: fig.savefig(a.out/f'Figure5_pairwise_ROC_24h.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
fig.savefig(a.out/'Figure5_pairwise_ROC_24h.tiff',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
fig.savefig(a.out/'JTIM_Figure_5.jpg',dpi=300,bbox_inches='tight')
plt.close(fig)

# ---------------- Figure S8 and S20 CSV (code from plot_windows.py)
ci=pd.read_csv(a.src/'auc_confidence_intervals.csv')
def ci_text(h,s,k):
 r=ci[(ci.hours==h)&(ci.split==s)&(ci.metric==k)].iloc[0]
 return f'{r.auc:.3f} ({r.lower95:.3f}–{r.upper95:.3f})'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'legend.fontsize':8,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,2,figsize=(7.2,5.8),sharex=True,sharey=True)
for ax,key,title,letter in zip(axs.flat,['auc_DR_PW','auc_ovo_macro','auc_DR_RR','auc_RR_PW'],
                             ['DR versus PW','Overall macro OvO','DR versus RR','RR versus PW'],'abcd'):
    for split,color,marker,label in [('internal','#37668C','o','Internal'),('external','#B6693E','s','External AUMC')]:
        z=ci[(ci.metric==key)&(ci.split==split)].sort_values('hours')
        ax.fill_between(z.hours.to_numpy(),z.lower95.to_numpy(),z.upper95.to_numpy(),color=color,alpha=.12,lw=0)
        ax.plot(z.hours,z.auc,color=color,marker=marker,ms=3.8,lw=1.5,label=label)
    ax.axvline(24,color='#888888',ls='--',lw=.8,zorder=0)
    ax.set_title(title,loc='left',pad=7);ax.text(-.15,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=11)
    ax.set_xticks([6,12,18,24,30,36]);ax.set_ylim(.5,1.0);ax.set_yticks(np.arange(.5,1.01,.1))
    ax.grid(axis='y',color='#eeeeee',lw=.5)
for ax in axs[:,0]:ax.set_ylabel('AUC')
for ax in axs[-1,:]:ax.set_xlabel('Hours after SA-AKI onset')
handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=2,bbox_to_anchor=(.54,1.01),frameon=False)
fig.subplots_adjust(left=.1,right=.98,bottom=.1,top=.91,wspace=.22,hspace=.36)
for ext in ['pdf','svg','png','tiff']:
    fig.savefig(a.out/f'window_prediction.{ext}',dpi=300 if ext=='png' else 600,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
plt.close(fig)
table=[]
for s in ['internal','external']:
 for h in [6,12,18,24,30,36]:
    table.append({'Cohort':s,'Window (h)':h,'Overall macro OvO':ci_text(h,s,'auc_ovo_macro'),
                  'DR vs RR':ci_text(h,s,'auc_DR_RR'),'DR vs PW':ci_text(h,s,'auc_DR_PW'),'RR vs PW':ci_text(h,s,'auc_RR_PW')})
pd.DataFrame(table).to_csv(a.out/'supplementary_table_auc_95CI.csv',index=False)
print('written', a.out)
