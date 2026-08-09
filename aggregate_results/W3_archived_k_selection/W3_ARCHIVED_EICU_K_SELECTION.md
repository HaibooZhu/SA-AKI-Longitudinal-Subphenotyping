# W3 归档 eICU 候选 K=2–5 动态复算

**状态：PASS_WITH_ARCHIVE_LIMITATIONS**

## 结论

- 直接读取冻结 `mixAK.RData` 中的 `mod2`–`mod5`，未使用 notebook 里手工抄写的结果。
- 归档只支持 K=2–5；没有 K=6–8 的模型对象。
- 按归档方法使用 deviance 与 |lag-1 autocorrelation|>0.85 的参数比例构造评分，K=3 排名第一。
- posterior median probability 保持原生 0–1 尺度，未再除以 2；这会修复概率图，但不改变 `which.max` 表型标签。
- K=4/5 仅有 500 个保留抽样，而 K=2/3 有 2,000 个，故高 K 比较证据较弱。
- 该归档仅能动态复算 eICU；另外两库缺少同等级 model object，不能声称三库均完成了相同的 K 诊断复算。

## 动态复算结果

 K patients retained_draws mean_deviance high_absolute_lag1_fraction
 2     1417           2000      556131.1                     0.00000
 3     1417           2000      555519.6                     0.00000
 4     1417            500      702017.7                     0.65625
 5     1417            500      702157.1                     0.50000
 uncertain_patients uncertain_fraction median_max_posterior_probability
                 44         0.03105152                        0.9992171
                206         0.14537756                        0.9643289
                 58         0.04093155                        0.9996607
                134         0.09456598                        0.9822235
        probability_scale scaled_deviance scaled_lag1_failure
 native 0-1 (no division)     0.004169909           0.0000000
 native 0-1 (no division)     0.000000000           0.0000000
 native 0-1 (no division)     0.999049349           1.0000000
 native 0-1 (no division)     1.000000000           0.7619048
 historical_selection_score rank_by_historical_score selected
                0.004169909                        2    FALSE
                0.000000000                        1     TRUE
                1.413541510                        4    FALSE
                1.257178932                        3    FALSE

冻结归档 SHA-256：`7fbd4265968a68b99c6948a86ef997e09c5e55daacaa4716ebd3af0dcaf22193`。
