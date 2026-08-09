#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_mixak_k3_sensitivity.R INPUT_CSV OUTPUT_DIR SCENARIO [SEED] [CONFIG_JSON]")
}

input_csv <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
scenario <- args[[3]]
seed <- if (length(args) >= 4) as.integer(args[[4]]) else 20260805L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
local_library <- normalizePath(
  file.path(dirname(script_path), "..", "00_config", "R_libs"),
  mustWork = FALSE
)
if (dir.exists(local_library)) {
  .libPaths(c(local_library, .libPaths()))
}

config_path <- if (length(args) >= 5) {
  normalizePath(args[[5]], mustWork = TRUE)
} else {
  normalizePath(
    file.path(dirname(script_path), "..", "00_config", "mixak_canonical_config.json"),
    mustWork = TRUE
  )
}

library(mixAK)
library(coda)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package jsonlite is required to read the canonical mixAK configuration")
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package digest is required to checksum the canonical mixAK configuration")
}
config <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
run_config <- config$revision_sensitivity
assignment_config <- config$posterior_assignment
diagnostic_config <- config$diagnostics

if (run_config$k != 3) {
  stop("This revision sensitivity runner requires canonical K=3")
}
if (run_config$chains != 1) {
  stop("This runner implements the one-chain archived reproduction profile")
}
if (assignment_config$probability_scaling != "none") {
  stop("Posterior probabilities must remain on their native 0-1 scale")
}

set.seed(seed)
df <- read.csv(input_csv, check.names = FALSE)
features <- as.character(run_config$features)
required <- c("stay_id", "time", features)
missing_columns <- setdiff(required, names(df))
if (length(missing_columns) > 0) {
  stop(sprintf("Missing columns: %s", paste(missing_columns, collapse = ", ")))
}
df <- df[complete.cases(df[, required]), required]
df <- df[order(df$stay_id, df$time), ]

model <- GLMM_MCMC(
  y = df[, features],
  dist = rep("gaussian", 1),
  id = df$stay_id,
  x = list(
    bun = "empty",
    creatinine = "empty",
    urineoutput = "empty",
    crea_divide_basecrea = "empty"
  ),
  z = list(
    bun = df$time,
    creatinine = df$time,
    urineoutput = df$time,
    crea_divide_basecrea = df$time
  ),
  random.intercept = rep(TRUE, length(features)),
  prior.b = list(Kmax = as.integer(run_config$k)),
  nMCMC = c(
    burn = as.integer(run_config$burn),
    keep = as.integer(run_config$keep),
    thin = as.integer(run_config$thin),
    info = as.integer(run_config$info)
  ),
  parallel = TRUE,
  PED = FALSE
)
model <- NMixRelabel(model, type = "stephens", keep.comp.prob = TRUE)
saveRDS(model, file.path(output_dir, paste0(scenario, "_model.rds")), compress = FALSE)

# quant.comp.prob is already on the 0-1 probability scale.
median_probabilities <- model$quant.comp.prob[["50%"]]
if (!identical(dim(median_probabilities), c(length(unique(df$stay_id)), as.integer(run_config$k)))) {
  stop("Posterior median probability matrix has unexpected dimensions")
}
if (any(!is.finite(median_probabilities)) || any(median_probabilities < 0) || any(median_probabilities > 1)) {
  stop("Posterior median probabilities fall outside the 0-1 scale")
}
group_median <- apply(median_probabilities, 1, which.max)
max_median_probability <- apply(median_probabilities, 1, max)
patient_ids <- unique(df$stay_id)
if (length(patient_ids) != length(group_median)) {
  stop(sprintf(
    "Patient/order mismatch: %d unique input IDs versus %d posterior rows",
    length(patient_ids), length(group_median)
  ))
}

hpd <- HPDinterval(mcmc(model$comp.prob), prob = assignment_config$hpd_probability)
probability_names <- rownames(hpd)
matches <- regexec("^P\\(([0-9]+),([0-9]+)\\)$", probability_names)
parsed <- regmatches(probability_names, matches)
if (any(lengths(parsed) != 3)) {
  stop("Could not parse posterior probability column names for HPD alignment")
}
indices <- do.call(rbind, lapply(parsed, function(value) as.integer(value[2:3])))
hpd_lower <- matrix(
  NA_real_,
  nrow = length(patient_ids),
  ncol = as.integer(run_config$k)
)
hpd_lower[cbind(indices[, 1], indices[, 2])] <- hpd[, "lower"]
if (anyNA(hpd_lower)) {
  stop("HPD lower-bound matrix is incomplete after index-based alignment")
}
group_hpd <- group_median
for (cluster in seq_len(as.integer(run_config$k))) {
  group_hpd[
    group_hpd == cluster &
      hpd_lower[, cluster] <= assignment_config$hpd_lower_threshold
  ] <- as.integer(assignment_config$uncertain_label)
}

assignments <- data.frame(
  stay_id = patient_ids,
  group_median = as.integer(group_median),
  group_hpd = as.integer(group_hpd),
  max_median_probability = as.numeric(max_median_probability),
  stringsAsFactors = FALSE
)
write.csv(assignments, file.path(output_dir, paste0(scenario, "_assignments.csv")), row.names = FALSE)

mu_chains <- NMixChainComp(model, relabel = TRUE, param = "mu_b")
lag_value <- as.integer(diagnostic_config$lag)
lag1 <- apply(
  mu_chains,
  2,
  function(x) as.numeric(autocorr(mcmc(x), lags = lag_value))
)
lag1_threshold <- as.numeric(diagnostic_config$absolute_autocorrelation_threshold)
high_lag1_fraction <- mean(abs(lag1) > lag1_threshold)
maximum_high_fraction <- as.numeric(diagnostic_config$maximum_fraction_above_threshold)
convergence_status <- if (high_lag1_fraction <= maximum_high_fraction) {
  "PASS_SINGLE_CHAIN_DIAGNOSTIC"
} else {
  "CONVERGENCE_WARNING"
}
diagnostics <- data.frame(
  scenario = scenario,
  seed = seed,
  config_sha256 = digest::digest(file = config_path, algo = "sha256"),
  k = as.integer(run_config$k),
  chains = as.integer(run_config$chains),
  burn = as.integer(run_config$burn),
  keep = as.integer(run_config$keep),
  thin = as.integer(run_config$thin),
  patients = length(unique(df$stay_id)),
  rows = nrow(df),
  mean_deviance = mean(model$Deviance),
  lag1_absolute_threshold = lag1_threshold,
  high_lag1_fraction = high_lag1_fraction,
  convergence_status = convergence_status,
  posterior_probability_min = min(median_probabilities),
  posterior_probability_max = max(median_probabilities),
  all_zero_posterior_median_rows = sum(rowSums(median_probabilities) == 0),
  uncertain_patients = sum(group_hpd == 4),
  uncertain_fraction = mean(group_hpd == 4),
  stringsAsFactors = FALSE
)
write.csv(diagnostics, file.path(output_dir, paste0(scenario, "_diagnostics.csv")), row.names = FALSE)

status <- list(
  scenario = scenario,
  status = convergence_status,
  config_path = config_path,
  config_sha256 = diagnostics$config_sha256,
  posterior_probability_scaling = assignment_config$probability_scaling,
  uncertainty_rule = assignment_config$uncertainty_rule,
  single_chain_limitation = diagnostic_config$single_chain_limitation
)
jsonlite::write_json(
  status,
  file.path(output_dir, paste0(scenario, "_run_status.json")),
  pretty = TRUE,
  auto_unbox = TRUE
)

print(diagnostics)
