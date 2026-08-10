#!/usr/bin/env bash

# Launch the prespecified revision-only deep mixAK grid in an isolated workspace.
# Frozen project files are read-only inputs; all fit objects are written below OUTPUT_ROOT.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_WORK_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK_ROOT="${JTIM_WORK_ROOT:-${DEFAULT_WORK_ROOT}}"
RSCRIPT_BIN="${JTIM_RSCRIPT_BIN:-${WORK_ROOT}/env/bin/Rscript}"
RUNNER="${WORK_ROOT}/code/03_cluster_robustness/run_mixak_refit.R"
CONFIG="${WORK_ROOT}/code/00_config/mixak_canonical_config.json"
INPUT_ROOT="${WORK_ROOT}/inputs"
OUTPUT_ROOT="${WORK_ROOT}/outputs"
MAX_JOBS="${JTIM_MAX_JOBS:-24}"
SEEDS=(20260805 20260806 20260807)

mkdir -p \
  "${OUTPUT_ROOT}/deep_k" \
  "${OUTPUT_ROOT}/deep_robustness" \
  "${OUTPUT_ROOT}/uo_multiseed" \
  "${OUTPUT_ROOT}/job_status"

run_one() {
  local job_id="$1"
  local input_csv="$2"
  local output_dir="$3"
  local cohort="$4"
  local scenario="$5"
  local k="$6"
  local seed="$7"
  local features="$8"
  local log_path="${OUTPUT_ROOT}/job_status/${job_id}.log"
  local exit_path="${OUTPUT_ROOT}/job_status/${job_id}.exit"

  if "${RSCRIPT_BIN}" "${RUNNER}" \
    "${input_csv}" "${output_dir}" "${cohort}" "${scenario}" \
    "${k}" "${seed}" "${features}" true "${CONFIG}" deep_refit \
    >"${log_path}" 2>&1; then
    printf '0\n' >"${exit_path}"
  else
    printf '%s\n' "$?" >"${exit_path}"
  fi
}

wait_for_capacity() {
  while (( $(jobs -pr | wc -l) >= MAX_JOBS )); do
    wait -n || true
  done
}

launch() {
  wait_for_capacity
  run_one "$@" &
}

for cohort in mimic eicu aumc; do
  if [[ "${cohort}" == "aumc" ]]; then
    features="creatinine,urineoutput,crea_divide_basecrea"
  else
    features="bun,creatinine,urineoutput,crea_divide_basecrea"
  fi
  for k in 2 3; do
    for seed in "${SEEDS[@]}"; do
      job_id="${cohort}__primary_deep__K${k}__seed${seed}"
      launch \
        "${job_id}" \
        "${INPUT_ROOT}/primary/${cohort}.csv" \
        "${OUTPUT_ROOT}/deep_k" \
        "${cohort}" "primary_deep" "${k}" "${seed}" "${features}"
    done
  done
done

ROBUSTNESS_JOBS=(
  "mimic|complete_30_window_followup|mimic__complete_30_window_followup.csv|bun,creatinine,urineoutput,crea_divide_basecrea"
  "mimic|limited_forward_fill_complete_rows|mimic__limited_forward_fill_complete_rows.csv|bun,creatinine,urineoutput,crea_divide_basecrea"
  "eicu|complete_30_window_followup|eicu__complete_30_window_followup.csv|bun,creatinine,urineoutput,crea_divide_basecrea"
  "eicu|limited_forward_fill_complete_rows|eicu__limited_forward_fill_complete_rows.csv|bun,creatinine,urineoutput,crea_divide_basecrea"
  "aumc|exclude_documented_rrt|aumc__exclude_documented_rrt.csv|creatinine,urineoutput,crea_divide_basecrea"
  "aumc|complete_30_window_followup|aumc__complete_30_window_followup.csv|creatinine,urineoutput,crea_divide_basecrea"
)

for spec in "${ROBUSTNESS_JOBS[@]}"; do
  IFS='|' read -r cohort scenario filename features <<<"${spec}"
  for seed in "${SEEDS[@]}"; do
    job_id="${cohort}__${scenario}_deep__K3__seed${seed}"
    launch \
      "${job_id}" \
      "${INPUT_ROOT}/robustness/${filename}" \
      "${OUTPUT_ROOT}/deep_robustness" \
      "${cohort}" "${scenario}_deep" 3 "${seed}" "${features}"
  done
done

for scenario in documented_windows high_coverage; do
  for seed in "${SEEDS[@]}"; do
    job_id="eicu__${scenario}_deep__K3__seed${seed}"
    launch \
      "${job_id}" \
      "${INPUT_ROOT}/uo/eicu_mixak_${scenario}.csv" \
      "${OUTPUT_ROOT}/uo_multiseed" \
      eicu "${scenario}_deep" 3 "${seed}" \
      "bun,creatinine,urineoutput,crea_divide_basecrea"
  done
done

wait

expected=42
observed=$(find "${OUTPUT_ROOT}/job_status" -name '*.exit' -type f | wc -l)
failed=$(awk '$1 != 0 {count += 1} END {print count + 0}' "${OUTPUT_ROOT}"/job_status/*.exit)
printf 'expected=%s\nobserved=%s\nfailed=%s\n' \
  "${expected}" "${observed}" "${failed}" >"${OUTPUT_ROOT}/ALL_DONE"

if (( observed != expected || failed != 0 )); then
  exit 1
fi
