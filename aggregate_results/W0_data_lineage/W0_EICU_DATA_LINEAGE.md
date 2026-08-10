# W0 eICU 数据血缘与队列一致性报告

**总体状态：PASS_WITH_QUARANTINED_LEGACY_ARTIFACTS**

本报告将 `03.eICU_SAKI_trajCluster/sk_survival.csv` 中的 1,417 名患者定义为本轮修订的权威 eICU 队列。患者标识只在内存中用于集合比较；所有落盘文件均只含聚合计数和不可逆 SHA-256 摘要。

## 结论先行

- 共核查 4 个最终分析输入，其中 4 个通过患者集合与标签一致性规则。
- 旧版生存时间汇总含 1,970 名患者；相对权威队列，额外 553 人、缺少 0 人，重叠患者中有 200 个表型标签差异和 0 个 28 天死亡结局差异。
- 旧版风险因素合并文件含 1,748 名患者；相对权威队列，额外 331 人、缺少 0 人，重叠患者中有 148 个表型标签差异和 0 个 28 天死亡结局差异。
- 分类器上游特征快照含 1,748 名患者；相对权威队列，额外 331 人、缺少 0 人，重叠患者中有 148 个表型标签差异和 0 个 28 天死亡结局差异。
- 分类器内部输入的队列归属检查发现 0 个跨数据库歧义 ID 和 0 个无法归属的内部拆分 ID；任一计数非零都会使审计失败。
- 因此，所有标记为 `NOT_AUTHORIZED_FOR_FINAL_RESULTS` 的旧派生文件只可用于追溯；最终分析必须从权威队列及经核验的上游变量重新派生。

## 工件清单

| artifact                               |   unique_patient_count | trust_decision                   | expected_relation   |
|:---------------------------------------|-----------------------:|:---------------------------------|:--------------------|
| eicu_authoritative_survival_1417       |                   1417 | REFERENCE                        | REFERENCE           |
| eicu_final_cluster_input               |                   1417 | PASS                             | MATCH               |
| eicu_final_imputed_matrix              |                   1417 | PASS                             | MATCH               |
| eicu_cross_cohort_longitudinal_source  |                   1417 | PASS                             | MATCH               |
| eicu_legacy_survival_time_summary      |                   1970 | NOT_AUTHORIZED_FOR_FINAL_RESULTS | TRACE_ONLY          |
| eicu_legacy_risk_factor_union          |                   1748 | NOT_AUTHORIZED_FOR_FINAL_RESULTS | TRACE_ONLY          |
| eicu_legacy_classifier_feature_source  |                   1748 | NOT_AUTHORIZED_FOR_FINAL_RESULTS | TRACE_ONLY          |
| eicu_classifier_actual_internal_inputs |                   1417 | PASS                             | MATCH               |

## 与权威队列的逐项比较

| artifact                               | relation                  |   overlap_n |   authoritative_only_n |   artifact_only_n |   group_mismatch_n |   mortality_mismatch_n |
|:---------------------------------------|:--------------------------|------------:|-----------------------:|------------------:|-------------------:|-----------------------:|
| eicu_authoritative_survival_1417       | MATCH                     |        1417 |                      0 |                 0 |                  0 |                      0 |
| eicu_final_cluster_input               | MATCH                     |        1417 |                      0 |                 0 |                  0 |                      0 |
| eicu_final_imputed_matrix              | MATCH                     |        1417 |                      0 |                 0 |                  0 |                      0 |
| eicu_cross_cohort_longitudinal_source  | MATCH                     |        1417 |                      0 |                 0 |                  0 |                      0 |
| eicu_legacy_survival_time_summary      | SUPERSET_OF_AUTHORITATIVE |        1417 |                      0 |               553 |                200 |                      0 |
| eicu_legacy_risk_factor_union          | SUPERSET_OF_AUTHORITATIVE |        1417 |                      0 |               331 |                148 |                      0 |
| eicu_legacy_classifier_feature_source  | SUPERSET_OF_AUTHORITATIVE |        1417 |                      0 |               331 |                148 |                      0 |
| eicu_classifier_actual_internal_inputs | MATCH                     |        1417 |                      0 |                 0 |                  0 |                      0 |

## 放行规则

1. 最终分析输入必须与权威患者集合完全一致。
2. 含表型标签的最终输入不得存在患者内冲突或与权威标签不一致。
3. 含 28 天死亡结局的最终结局输入不得存在患者内冲突或与权威结局不一致。
4. 标记为 `NOT_AUTHORIZED_FOR_FINAL_RESULTS` 的旧快照只可用于追溯问题，不得直接进入修订分析。

详细文件摘要见 `eicu_artifact_manifest.csv`，逐项集合比较见 `eicu_pairwise_lineage.csv`，机器可读状态见 `eicu_lineage_status.json`。
