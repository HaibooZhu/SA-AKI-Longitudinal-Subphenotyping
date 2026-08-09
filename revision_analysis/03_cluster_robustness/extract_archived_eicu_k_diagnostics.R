#!/usr/bin/env Rscript

# Recompute aggregate K=2-5 diagnostics from the frozen eICU mixAK archive.
# Frozen model objects are read only; no row-level assignments are exported.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
if (length(args) < 1) {
  stop("Usage: extract_archived_eicu_k_diagnostics.R ARCHIVE_RDATA [OUTPUT_DIR]")
}
archive_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path("results", "revision", "W3_archived_k_selection")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

local_library <- normalizePath(
  file.path(dirname(script_path), "..", "00_config", "R_libs"),
  mustWork = FALSE
)
if (dir.exists(local_library)) {
  .libPaths(c(local_library, .libPaths()))
}

suppressPackageStartupMessages(library(mixAK))
suppressPackageStartupMessages(library(coda))

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package digest is required")
  }
  digest::digest(file = path, algo = "sha256")
}

load(archive_path)
required_objects <- paste0("mod", 2:5)
missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects) > 0) {
  stop(sprintf("Archive lacks required objects: %s", paste(missing_objects, collapse = ", ")))
}

diagnose_model <- function(model_name) {
  model <- get(model_name, envir = .GlobalEnv)
  k <- as.integer(model$prior.b$Kmax)
  median_probabilities <- model$quant.comp.prob[["50%"]]
  if (!identical(ncol(median_probabilities), k)) {
    stop(sprintf("%s: posterior probability column count does not match K", model_name))
  }
  if (any(!is.finite(median_probabilities)) ||
      any(median_probabilities < 0) ||
      any(median_probabilities > 1)) {
    stop(sprintf("%s: posterior medians fall outside native 0-1 scale", model_name))
  }
  group_median <- apply(median_probabilities, 1, which.max)

  hpd <- HPDinterval(mcmc(model$comp.prob), prob = 0.95)
  probability_names <- rownames(hpd)
  matches <- regexec("^P\\(([0-9]+),([0-9]+)\\)$", probability_names)
  parsed <- regmatches(probability_names, matches)
  if (any(lengths(parsed) != 3)) {
    stop(sprintf("%s: could not parse posterior-probability indices", model_name))
  }
  indices <- do.call(rbind, lapply(parsed, function(value) as.integer(value[2:3])))
  hpd_lower <- matrix(NA_real_, nrow = nrow(median_probabilities), ncol = k)
  hpd_lower[cbind(indices[, 1], indices[, 2])] <- hpd[, "lower"]
  if (anyNA(hpd_lower)) {
    stop(sprintf("%s: incomplete HPD lower-bound matrix", model_name))
  }
  uncertain <- rep(FALSE, length(group_median))
  for (cluster in seq_len(k)) {
    uncertain[group_median == cluster & hpd_lower[, cluster] <= 0.5] <- TRUE
  }

  mu_chains <- NMixChainComp(model, relabel = TRUE, param = "mu_b")
  lag1 <- apply(
    mu_chains,
    2,
    function(x) as.numeric(autocorr(mcmc(x), lags = 1))
  )
  data.frame(
    K = k,
    patients = nrow(median_probabilities),
    retained_draws = as.integer(model$nMCMC[[2]]),
    mean_deviance = mean(model$Deviance),
    high_absolute_lag1_fraction = mean(abs(lag1) > 0.85),
    uncertain_patients = sum(uncertain),
    uncertain_fraction = mean(uncertain),
    median_max_posterior_probability = median(apply(median_probabilities, 1, max)),
    probability_scale = "native 0-1 (no division)",
    stringsAsFactors = FALSE
  )
}

diagnostics <- do.call(rbind, lapply(required_objects, diagnose_model))
diagnostics <- diagnostics[order(diagnostics$K), ]

scale_zero_one <- function(values) {
  if (max(values) == min(values)) return(rep(0, length(values)))
  (values - min(values)) / (max(values) - min(values))
}
diagnostics$scaled_deviance <- scale_zero_one(diagnostics$mean_deviance)
diagnostics$scaled_lag1_failure <- scale_zero_one(diagnostics$high_absolute_lag1_fraction)
diagnostics$historical_selection_score <- sqrt(
  diagnostics$scaled_deviance^2 + diagnostics$scaled_lag1_failure^2
)
diagnostics$rank_by_historical_score <- rank(
  diagnostics$historical_selection_score,
  ties.method = "min"
)
diagnostics$selected <- diagnostics$rank_by_historical_score == 1

selected_k <- diagnostics$K[diagnostics$selected]
if (!identical(selected_k, 3L)) {
  stop(sprintf("Archived diagnostics selected unexpected K: %s", paste(selected_k, collapse = ",")))
}
write.csv(diagnostics, file.path(output_dir, "archived_eicu_k2_k5_diagnostics.csv"), row.names = FALSE)

png(file.path(output_dir, "Figure_S2_archived_eicu_k_selection.png"), width = 1800, height = 1050, res = 180)
par(mar = c(5, 5, 3, 5) + 0.1, family = "sans")
plot(
  diagnostics$K,
  diagnostics$historical_selection_score,
  type = "b", pch = 19, lwd = 2.5,
  col = "#2F6B9A", xaxt = "n",
  xlab = "Candidate number of clusters (K)",
  ylab = "Archived selection score (lower is better)",
  ylim = c(0, max(diagnostics$historical_selection_score) * 1.12)
)
axis(1, at = diagnostics$K)
grid(nx = NA, ny = NULL, col = "#D9E2EA", lty = 1)
points(3, diagnostics$historical_selection_score[diagnostics$K == 3], pch = 21, cex = 1.7, bg = "#E67E22", col = "#A84F0A", lwd = 2)
text(3, diagnostics$historical_selection_score[diagnostics$K == 3], labels = "Selected K=3", pos = 3, cex = 1.0)
par(new = TRUE)
plot(
  diagnostics$K,
  diagnostics$uncertain_fraction,
  type = "b", pch = 17, lwd = 2, lty = 2,
  col = "#6A8E3A", axes = FALSE, xlab = "", ylab = "",
  ylim = c(0, max(diagnostics$uncertain_fraction) * 1.25)
)
axis(4, col.axis = "#4E6F27")
mtext("Uncertain assignment fraction", side = 4, line = 3, col = "#4E6F27")
legend(
  "topleft",
  legend = c("Selection score", "Uncertain fraction"),
  col = c("#2F6B9A", "#6A8E3A"),
  pch = c(19, 17), lty = c(1, 2), lwd = 2,
  bty = "n"
)
title(main = "Archived eICU K=2-5 model diagnostics")
dev.off()

pdf(file.path(output_dir, "Figure_S2_archived_eicu_k_selection.pdf"), width = 10, height = 6)
par(mar = c(5, 5, 3, 5) + 0.1, family = "sans")
plot(
  diagnostics$K,
  diagnostics$historical_selection_score,
  type = "b", pch = 19, lwd = 2.5,
  col = "#2F6B9A", xaxt = "n",
  xlab = "Candidate number of clusters (K)",
  ylab = "Archived selection score (lower is better)",
  ylim = c(0, max(diagnostics$historical_selection_score) * 1.12)
)
axis(1, at = diagnostics$K)
grid(nx = NA, ny = NULL, col = "#D9E2EA", lty = 1)
points(3, diagnostics$historical_selection_score[diagnostics$K == 3], pch = 21, cex = 1.7, bg = "#E67E22", col = "#A84F0A", lwd = 2)
text(3, diagnostics$historical_selection_score[diagnostics$K == 3], labels = "Selected K=3", pos = 3, cex = 1.0)
par(new = TRUE)
plot(
  diagnostics$K,
  diagnostics$uncertain_fraction,
  type = "b", pch = 17, lwd = 2, lty = 2,
  col = "#6A8E3A", axes = FALSE, xlab = "", ylab = "",
  ylim = c(0, max(diagnostics$uncertain_fraction) * 1.25)
)
axis(4, col.axis = "#4E6F27")
mtext("Uncertain assignment fraction", side = 4, line = 3, col = "#4E6F27")
legend(
  "topleft",
  legend = c("Selection score", "Uncertain fraction"),
  col = c("#2F6B9A", "#6A8E3A"),
  pch = c(19, 17), lty = c(1, 2), lwd = 2,
  bty = "n"
)
title(main = "Archived eICU K=2-5 model diagnostics")
dev.off()

report_lines <- c(
  "# W3 归档 eICU 候选 K=2–5 动态复算",
  "",
  "**状态：PASS_WITH_ARCHIVE_LIMITATIONS**",
  "",
  "## 结论",
  "",
  "- 直接读取冻结 `mixAK.RData` 中的 `mod2`–`mod5`，未使用 notebook 里手工抄写的结果。",
  "- 归档只支持 K=2–5；没有 K=6–8 的模型对象。",
  "- 按归档方法使用 deviance 与 |lag-1 autocorrelation|>0.85 的参数比例构造评分，K=3 排名第一。",
  "- posterior median probability 保持原生 0–1 尺度，未再除以 2；这会修复概率图，但不改变 `which.max` 表型标签。",
  "- K=4/5 仅有 500 个保留抽样，而 K=2/3 有 2,000 个，故高 K 比较证据较弱。",
  "- 该归档仅能动态复算 eICU；另外两库缺少同等级 model object，不能声称三库均完成了相同的 K 诊断复算。",
  "",
  "## 动态复算结果",
  "",
  paste(capture.output(print(diagnostics, row.names = FALSE)), collapse = "\n"),
  "",
  sprintf("冻结归档 SHA-256：`%s`。", sha256_file(archive_path)),
  ""
)
writeLines(report_lines, file.path(output_dir, "W3_ARCHIVED_EICU_K_SELECTION.md"), useBytes = TRUE)

status <- list(
  overall_status = "PASS_WITH_ARCHIVE_LIMITATIONS",
  archived_candidate_k = 2:5,
  selected_k = 3,
  unsupported_k = 6:8,
  probability_scaling = "none",
  patient_level_output_written = FALSE,
  archive_sha256 = sha256_file(archive_path)
)
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package jsonlite is required")
}
jsonlite::write_json(
  status,
  file.path(output_dir, "archived_eicu_k_selection_status.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
print(diagnostics)
