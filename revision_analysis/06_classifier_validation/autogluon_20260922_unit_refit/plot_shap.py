"""Plot only aggregate TreeSHAP values verified with AutoGluon CSR inputs."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parents[3]/'02_revision_outputs/reports/20260922_unit_harmonized_refit'
x=pd.read_csv(OUT/'shap_mean_abs_by_feature_24h.csv').sort_values('total',ascending=False).head(15)
assert np.allclose(x[['DR','RR','PW']].sum(axis=1),x.total)
def label(f):
 names={'crea_divide_basecrea':'SCr/bSCr ratio','creatinine':'Serum creatinine','hemoglobin':'Hemoglobin','urineoutput':'Urine output','resp_rate':'Respiratory rate','calcium':'Calcium','bilirubin':'Bilirubin'}
 for base,name in names.items():
  if f.startswith(base+'_'):
   suffix=f[len(base)+1:]
   if suffix.startswith('diff_'):
    a,b=suffix.split('_')[1:];return f'{name}, change from window {a} to {b}'
   return name+', '+{'max':'maximum','min':'minimum','mean':'mean'}[suffix]
 raise ValueError(f)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(8.4,6.1)); y=np.arange(len(x)); left=np.zeros(len(x))
for name,color in [('DR','#1f77b4'),('RR','#2ca02c'),('PW','#d62728')]:
 ax.barh(y,x[name],left=left,height=.65,color=color,label=name);left+=x[name].to_numpy()
ax.errorbar(x.total,y,xerr=x.fold_sd_total,fmt='none',ecolor='black',elinewidth=.8,capsize=2)
for yy,total,sd in zip(y,x.total,x.fold_sd_total):ax.text(total+sd+.011,yy,f'{total:.2f}',va='center',fontsize=8)
ax.set_yticks(y,[label(f) for f in x.feature]);ax.invert_yaxis();ax.set_xlim(0,max(x.total+x.fold_sd_total)*1.17)
ax.set_xlabel('Mean |SHAP value| (raw-margin scale), averaged across the five fold models',fontsize=8)
ax.legend(loc='lower right',frameon=False);fig.tight_layout()
for ext in ['pdf','svg','png','tiff']:
 fig.savefig(OUT/f'FigureS7_SHAP_24h.{ext}',dpi=300 if ext=='png' else 600,bbox_inches='tight',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
