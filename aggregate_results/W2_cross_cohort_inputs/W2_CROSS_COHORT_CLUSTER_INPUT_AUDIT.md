# W2 跨队列聚类输入、单位与基线肌酐审计

**状态：PASS_WITH_DISCLOSED_STRUCTURAL_HETEROGENEITY**

## 审计结论

- 三个冻结聚类矩阵均无重复 patient-time 行，每位患者只有一个归档表型标签。
- MIMIC-IV 与 eICU-CRD 的可执行聚类输入包含 BUN、肌酐、6 小时尿量和肌酐/基线肌酐比值 4 个肾脏变量。
- AUMC 的归档矩阵不含 BUN，实际使用后 3 个变量。因此修订稿必须描述为“两个四变量模型和一个三变量外部验证模型”，不能声称三个数据库使用完全相同的四变量输入。
- 尿量字段的可审计含义为每个 6 小时窗口内的尿量体积（mL）；不是 24 小时总量，也没有按体重归一化。
- 三队列均使用单一归档基线肌酐字段，入模前要求 0.5–<1.5 mg/dL；缺失或超出范围者被排除，而不是插补。AUMC 原始单位为 µmol/L，按 0.01131 转换为 mg/dL。
- 冻结聚类归档中没有第二套可独立重建的基线肌酐定义。因此本次修订可澄清并审计既有定义，但不能伪造“替代基线定义”敏感性；这项限制应在回复和正文中明确承认。
- 肌酐/基线肌酐比值可由归档肌酐与基线值按两位小数重建，重建统计见下表。

## 矩阵结构

| cohort   | matrix_sha256                                                    |   rows |   patients |   time_min |   time_max |   unique_time_values |   patient_time_duplicates |   minimum_rows_per_patient |   median_rows_per_patient |   maximum_rows_per_patient |   renal_feature_count | renal_features                                  | BUN_in_executable_matrix   |
|:---------|:-----------------------------------------------------------------|-------:|-----------:|-----------:|-----------:|---------------------:|--------------------------:|---------------------------:|--------------------------:|---------------------------:|----------------------:|:------------------------------------------------|:---------------------------|
| MIMIC-IV | aca43c03eda33ae09354f474a10ecbfd3084351cfb60497b484793fb6ef4314c |  85021 |       4713 |         -2 |         28 |                   30 |                         0 |                          5 |                        16 |                         30 |                     4 | bun;creatinine;urineoutput;crea_divide_basecrea | True                       |
| eICU-CRD | a79251d11608c48c8c02aaadfe0c16c4a54dbf6ca5efaebd5efa34b550bface4 |  29246 |       1417 |         -2 |         28 |                   30 |                         0 |                          6 |                        21 |                         30 |                     4 | bun;creatinine;urineoutput;crea_divide_basecrea | True                       |
| AUMC     | 7bbc9080a59cda324ed14cd7126871a7c165b34fc0660bc4bd3dcb7508af96af |  48517 |       2183 |         -2 |         28 |                   30 |                         0 |                          5 |                        29 |                         30 |                     3 | creatinine;urineoutput;crea_divide_basecrea     | False                      |

## 实际进入聚类的肾脏变量

| feature              | AUMC   | MIMIC-IV   | eICU-CRD   |
|:---------------------|:-------|:-----------|:-----------|
| bun                  | False  | True       | True       |
| crea_divide_basecrea | True   | True       | True       |
| creatinine           | True   | True       | True       |
| urineoutput          | True   | True       | True       |

## 基线肌酐与比值重建

| cohort   | baseline_source_unit   |   conversion_to_mg_dl |   matrix_patients |   matrix_patients_linked_to_in_range_baseline |   ratio_exact_after_rounding_fraction |   ratio_absolute_error_max |
|:---------|:-----------------------|----------------------:|------------------:|----------------------------------------------:|--------------------------------------:|---------------------------:|
| MIMIC-IV | mg/dL                  |               1       |              4713 |                                          4713 |                                     1 |                          0 |
| eICU-CRD | mg/dL                  |               1       |              1417 |                                          1417 |                                     1 |                          0 |
| AUMC     | micromol/L             |               0.01131 |              2183 |                                          2183 |                                     1 |                          0 |

## 修订口径

1. Methods 按队列分别列出实际变量，而不是笼统写“所有数据库均为四变量”。
2. 将基线肌酐缺失/范围排除写入纳排流程，并在局限性中说明这限制了对既往肾功能异常人群的外推。
3. 将 AUMC 缺少 BUN 作为跨库结构差异，不把它包装成完全相同模型的独立复制。
4. 单位和数值范围汇总见 `cluster_feature_distribution_audit.csv`；全部结果仅为聚合统计，不含患者标识。
