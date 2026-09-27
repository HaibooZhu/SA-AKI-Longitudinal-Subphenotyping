"""Build an aggregate-only comparison with the archived primary classifier."""
import json
from pathlib import Path
import pandas as pd

repo=Path(__file__).resolve().parents[3]
out=repo/'02_revision_outputs/reports/20260921_autogluon_original_labels'
old=repo/'02_revision_outputs/reports/W4_classifier_validation'
new=pd.read_csv(out/'test_metrics.csv')
prior=pd.read_csv(old/'primary_and_comparator_metrics.csv')
prior=prior[prior.model.eq('archived_AutoGluon_XGBoost_BAG_L2')]
pairs=pd.read_csv(old/'primary_and_comparator_pairwise_auc.csv')
pairs=pairs[pairs.model.eq('archived_AutoGluon_XGBoost_BAG_L2')]
rows=[]
for r in prior.to_dict('records'):
    s='external' if r['cohort']=='external_AUMC' else 'internal'
    row={'version':'原论文保存模型','role':'original','split':s,'model':r['model'],
         **{k:r[k] for k in ['n','auc_ovo_macro','accuracy','balanced_accuracy','f1_macro','log_loss']}}
    row['brier_multiclass']=r['multiclass_brier']
    for a,b in [('DR vs RR','DR_RR'),('DR vs PW','DR_PW'),('RR vs PW','RR_PW')]:
        row['auc_'+b]=pairs[(pairs.cohort==r['cohort'])&(pairs.comparison==a)].iloc[0].auc_pair_normalized
    rows.append(row)
for r in new.to_dict('records'):
    rows.append({'version':'本轮 AutoGluon 1.6.3',**r})
comparison=pd.DataFrame(rows)
comparison.to_csv(out/'original_vs_new_metrics.csv',index=False)
lb=pd.read_csv(out/'training_cv_leaderboard.csv')
verification=json.loads((out/'verification.json').read_text())
locked=json.loads((out/'locked_selection.json').read_text())['selected']
lines=['# 原论文标签下的 AutoGluon 调参结果','',
       '本轮已完成服务器训练、测试集评价和结果复核。正文、补充材料、回复信及原始分析文件均未修改。','',
       '## 测试结果','',
       '总体 AUC 为 macro one-vs-one；成对 AUC 使用成对概率归一化。模型先由训练集五折结果选定，再评价测试集。','',
       '| 模型 | 测试集 | 总体 AUC | DR–RR | DR–PW | RR–PW | 准确率 |',
       '|---|---|---:|---:|---:|---:|---:|']
for r in rows:
    label={'original':'原论文模型','best_XGBoost':'本轮最佳 XGBoost','best_overall':'本轮总体最佳'}[r['role']]
    split={'internal':'内部','external':'AUMC 外部'}[r['split']]
    lines.append(f"| {label} | {split} | {r['auc_ovo_macro']:.3f} | {r['auc_DR_RR']:.3f} | {r['auc_DR_PW']:.3f} | {r['auc_RR_PW']:.3f} | {r['accuracy']:.3f} |")
lines += ['', '## 选择和实际配置','',f"- 最佳 XGBoost：`{locked['best_XGBoost']}`。",f"- 总体最佳：`{locked['best_overall']}`。",
          f'- 共完成 {len(verification["models"])} 个基础参数配置，所有基础模型均核实为 5 折、1 次重复。',
          '- 保留论文原分组；训练 4,903 人、内部测试 1,227 人、AUMC 外部测试 2,183 人。',
          '- 仅使用 SA-AKI 后 0–24 小时特征；共 27 个变量、243 个衍生特征。',
          '- 本轮未先做全训练集的结局相关特征筛选；243 个特征均进入 AutoGluon。与原模型相比，软件版本、特征处理和拟合方案均有变化，不能把全部差异归因于升级版本或某一个参数。',
          '- 搜索种子 20260921；五折袋装、无堆叠、允许加权集成；16 CPU、无 GPU；独立 Conda 环境，原环境保留。',
          '- AUMC 未参与本轮参数或模型选择。本轮测试结果不能用于继续筛选“外部最好”的参数后再声称独立验证。',
          '', '## 复核与文件','',
          '- 三个主分组文件哈希与此前审计一致；逐人标签与原分类器输入一致；三组患者互不重叠。',
          '- 所有特征与新实验包原构建脚本逐值比较通过，最大浮点差异 5.68e-14。',
          '- 已独立重算选定模型的交叉验证 AUC 和全部测试指标，核对五折实际完成、概率列顺序及原标签对应。',
          '- 患者级预测、模型和训练数据仅保存在服务器私有目录 `~/jtim_autogluon_20260921/`。',
          '- 本目录为聚合结果：`original_vs_new_metrics.csv`、`training_cv_leaderboard.csv`、`model_hyperparameters.json`、`verification.json`；来源和环境见 `input_manifest.json`、`environment.lock.txt`。',
          '', '## 与修稿的关系','',
          '这轮用于判断是否值得更新分类器结果，不自动替换论文。若决定采用，需按最终选定的模型同步方法、性能表和相关图件；无需因为本次实验重做聚类或 Table 1。',
          '', '版本与接口依据：[AutoGluon 1.6.3](https://pypi.org/project/autogluon.tabular/1.6.3/)；[官方训练接口](https://auto.gluon.ai/stable/api/autogluon.tabular.TabularPredictor.fit.html)。']
(out/'00_结果说明.md').write_text('\n'.join(lines)+'\n')
print(comparison[['role','split','auc_ovo_macro','auc_DR_PW','accuracy']].to_string(index=False))
