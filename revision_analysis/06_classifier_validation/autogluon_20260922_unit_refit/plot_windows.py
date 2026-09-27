from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parents[3]/'02_revision_outputs/reports/20260922_unit_harmonized_refit'
ci=pd.read_csv(OUT/'auc_confidence_intervals.csv')
def ci_text(h,s,k):
 r=ci[(ci.hours==h)&(ci.split==s)&(ci.metric==k)].iloc[0]
 return f'{r.auc:.3f} ({r.lower95:.3f}–{r.upper95:.3f})'
# Figure contract: quantitative grid; compare all four discrimination metrics,
# with DR-PW first. Fixed populations, full six windows, paired-bootstrap source.
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
    fig.savefig(OUT/f'window_prediction.{ext}',dpi=300 if ext=='png' else 600,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
plt.close(fig)

table=[]
for s in ['internal','external']:
 for h in [6,12,18,24,30,36]:
    table.append({'Cohort':s,'Window (h)':h,'Overall macro OvO':ci_text(h,s,'auc_ovo_macro'),
                  'DR vs RR':ci_text(h,s,'auc_DR_RR'),'DR vs PW':ci_text(h,s,'auc_DR_PW'),'RR vs PW':ci_text(h,s,'auc_RR_PW')})
pd.DataFrame(table).to_csv(OUT/'supplementary_table_auc_95CI.csv',index=False)
