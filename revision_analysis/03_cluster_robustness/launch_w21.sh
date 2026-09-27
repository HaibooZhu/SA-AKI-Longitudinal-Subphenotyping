#!/usr/bin/env bash
# W21 burn-in remediation launcher (remote).
#
# Single-factor re-estimation: identical data, patients, windows, responses,
# seeds, retained draws and thinning as the archived burn=50 deep grid; only the
# burn-in changes, from 50 to 20,000 iterations. Uses the same conda R 4.4.3 +
# mixAK 5.8 environment that produced the burn=50 grid, so the two grids differ
# in exactly one setting.
set -uo pipefail

ROOT="${JTIM_W21_ROOT:-$HOME/jtim_w21_burnin_20260823}"
SRC="${JTIM_SRC_ROOT:-$HOME/jtim_revision_20260809}"
RSCRIPT="${SRC}/env/bin/Rscript"
RUNNER="${ROOT}/code/run_mixak_converged_refit.R"
CONFIG="${ROOT}/code/mixak_converged_config.json"
OUT="${ROOT}/outputs"
LOG="${ROOT}/logs"
MAX_JOBS="${JTIM_MAX_JOBS:-42}"

mkdir -p "${OUT}" "${LOG}/job_status"

F4="bun,creatinine,urineoutput,crea_divide_basecrea"
F3="creatinine,urineoutput,crea_divide_basecrea"
SEEDS=(20260805 20260806 20260807)

build_joblist() {
  local jl="${ROOT}/joblist.txt"
  : > "${jl}"
  # Phase 1: primary K=2 and K=3, three cohorts, three initializations.
  for seed in "${SEEDS[@]}"; do
    for k in 2 3; do
      echo "mimic|primary_converged|${SRC}/inputs/primary/mimic.csv|${F4}|${k}|${seed}" >> "${jl}"
      echo "eicu|primary_converged|${SRC}/inputs/primary/eicu.csv|${F4}|${k}|${seed}"   >> "${jl}"
      echo "aumc|primary_converged|${SRC}/inputs/primary/aumc.csv|${F3}|${k}|${seed}"   >> "${jl}"
    done
  done
  # Phase 2: every screening-depth CAUTION scenario and both urine-output schemes.
  for seed in "${SEEDS[@]}"; do
    echo "eicu|documented_windows_converged|${SRC}/inputs/uo/eicu_mixak_documented_windows.csv|${F4}|3|${seed}" >> "${jl}"
    echo "eicu|high_coverage_converged|${SRC}/inputs/uo/eicu_mixak_high_coverage.csv|${F4}|3|${seed}" >> "${jl}"
    echo "eicu|complete_30_window_followup_converged|${SRC}/inputs/robustness/eicu__complete_30_window_followup.csv|${F4}|3|${seed}" >> "${jl}"
    echo "eicu|limited_forward_fill_complete_rows_converged|${SRC}/inputs/robustness/eicu__limited_forward_fill_complete_rows.csv|${F4}|3|${seed}" >> "${jl}"
    echo "aumc|exclude_documented_rrt_converged|${SRC}/inputs/robustness/aumc__exclude_documented_rrt.csv|${F3}|3|${seed}" >> "${jl}"
    echo "aumc|complete_30_window_followup_converged|${SRC}/inputs/robustness/aumc__complete_30_window_followup.csv|${F3}|3|${seed}" >> "${jl}"
    echo "mimic|complete_30_window_followup_converged|${SRC}/inputs/robustness/mimic__complete_30_window_followup.csv|${F4}|3|${seed}" >> "${jl}"
    echo "mimic|limited_forward_fill_complete_rows_converged|${SRC}/inputs/robustness/mimic__limited_forward_fill_complete_rows.csv|${F4}|3|${seed}" >> "${jl}"
  done
  echo "${jl}"
}

run_one() {
  local cohort="$1" scenario="$2" csv="$3" feats="$4" k="$5" seed="$6"
  local id="${cohort}__${scenario}__K${k}__seed${seed}"
  if OMP_NUM_THREADS="${JTIM_THREADS:-4}" OPENBLAS_NUM_THREADS="${JTIM_THREADS:-4}" \
      MKL_NUM_THREADS="${JTIM_THREADS:-4}" \
      "${RSCRIPT}" "${RUNNER}" "${csv}" "${OUT}" "${cohort}" "${scenario}" \
      "${k}" "${seed}" "${feats}" true "${CONFIG}" converged_refit \
      > "${LOG}/job_status/${id}.log" 2>&1; then
    printf '0\n' > "${LOG}/job_status/${id}.exit"
  else
    printf '%s\n' "$?" > "${LOG}/job_status/${id}.exit"
  fi
}

JOBLIST="$(build_joblist)"
TOTAL=$(grep -cve '^[[:space:]]*$' "${JOBLIST}")
printf 'started=%s\ntotal=%s\nmax_jobs=%s\n' "$(date -u +%FT%TZ)" "${TOTAL}" "${MAX_JOBS}" > "${LOG}/RUN_INFO"

while IFS='|' read -r cohort scenario csv feats k seed; do
  [[ -z "${cohort:-}" ]] && continue
  while (( $(jobs -pr | wc -l) >= MAX_JOBS )); do sleep 5; done
  run_one "${cohort}" "${scenario}" "${csv}" "${feats}" "${k}" "${seed}" &
done < "${JOBLIST}"
wait

observed=$(find "${LOG}/job_status" -name '*.exit' -type f | wc -l)
failed=$(awk '$1 != 0 {n += 1} END {print n + 0}' "${LOG}"/job_status/*.exit)
printf 'expected=%s\nobserved=%s\nfailed=%s\nfinished=%s\n' \
  "${TOTAL}" "${observed}" "${failed}" "$(date -u +%FT%TZ)" > "${LOG}/ALL_DONE"
cat "${LOG}/ALL_DONE"
