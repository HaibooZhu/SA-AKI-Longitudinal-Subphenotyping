# W3 mixAK 候选 K 与概率尺度溯源

**状态：PASS_WITH_REQUIRED_TABLE_CORRECTION**

## 结论

- 冻结的可执行 notebook、导出脚本与模型归档只支持候选 K=2, 3, 4, 5。
- 未发现任何 K=6、7、8 的可执行拟合路径、模型对象或日志，因此修订版 Table S1/Figure S2 不得继续保留 K=6–8。
- 冻结历史代码中发现 22 处 posterior median probability `/2` 写法。该常数缩放不改变 `which.max` 标签，但会破坏概率尺度及诊断图解释。
- 冻结历史/公开代码保持只读；修订 runner 使用原始 0–1 概率，并把配置文件 SHA-256 写入每次运行诊断。

## 允许进入修订稿的配置

- 历史可追溯比较：K=2–5。
- 当前尿量敏感性：固定 K=3，用于标签稳定性，不作为独立的 cluster-number selection 证据。
- 历史 K=4/5 仅保留 500 次迭代，弱于 K=2/3 的 2,000 次；这一限制必须如实披露。

机器可读证据见 `mixak_code_provenance.csv`、`mixak_rdata_provenance.csv` 与 `mixak_provenance_status.json`。
