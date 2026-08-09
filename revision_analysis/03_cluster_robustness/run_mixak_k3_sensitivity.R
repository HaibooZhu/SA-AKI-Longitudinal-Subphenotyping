#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_mixak_k3_sensitivity.R INPUT_CSV OUTPUT_DIR SCENARIO [SEED]")
}

input_csv <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
scenario <- args[[3]]
seed <- if (length(args) >= 4) as.integer(args[[4]]) else 20260805L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(mixAK)
library(coda)

set.seed(seed)
df <- read.csv(input_csv, check.names = FALSE)
features <- c("bun", "creatinine", "urineoutput", "crea_divide_basecrea")
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
  prior.b = list(Kmax = 3),
  nMCMC = c(burn = 50, keep = 2000, thin = 50, info = 50),
  parallel = TRUE,
  PED = FALSE
)
model <- NMixRelabel(model, type = "stephens", keep.comp.prob = TRUE)
saveRDS(model, file.path(output_dir, paste0(scenario, "_model.rds")), compress = FALSE)

# quant.comp.prob is already on the 0-1 probability scale.
median_probabilities <- model$quant.comp.prob[["50%"]]
group_median <- apply(median_probabilities, 1, which.max)
max_median_probability <- apply(median_probabilities, 1, max)
patient_ids <- unique(df$stay_id)
if (length(patient_ids) != length(group_median)) {
  stop(sprintf(
    "Patient/order mismatch: %d unique input IDs versus %d posterior rows",
    length(patient_ids), length(group_median)
  ))
}

hpd <- HPDinterval(mcmc(model$comp.prob))
hpd_lower <- matrix(hpd[, "lower"], ncol = 3, byrow = TRUE)
group_hpd <- group_median
for (cluster in 1:3) {
  group_hpd[group_hpd == cluster & hpd_lower[, cluster] <= 0.5] <- 4
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
lag1 <- apply(mu_chains, 2, function(x) autocorr(mcmc(x), lags = 1))
diagnostics <- data.frame(
  scenario = scenario,
  seed = seed,
  patients = length(unique(df$stay_id)),
  rows = nrow(df),
  mean_deviance = mean(model$Deviance),
  high_lag1_fraction = mean(lag1 > 0.85),
  uncertain_patients = sum(group_hpd == 4),
  uncertain_fraction = mean(group_hpd == 4),
  stringsAsFactors = FALSE
)
write.csv(diagnostics, file.path(output_dir, paste0(scenario, "_diagnostics.csv")), row.names = FALSE)

print(diagnostics)
