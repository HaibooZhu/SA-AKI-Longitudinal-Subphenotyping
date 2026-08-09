#!/usr/bin/env python3
"""Replay the archived AutoGluon 0.7.0 classifier and compare a simple model.

Run this script in the archived software environment (AutoGluon 0.7.0,
scikit-learn 1.2.2, XGBoost 1.7.4).  It writes aggregate metrics and calibration
bins only; it never exports patient identifiers or patient-level predictions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
import statsmodels.api as sm
import xgboost
from scipy.special import logit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    log_loss,
    precision_recall_fscore_support,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SEED = 20260805
CLASSES = np.array([1, 2, 3], dtype=int)
PHENOTYPE = {1: "DR", 2: "RR", 3: "PW"}
EXPECTED_VERSIONS = {
    "autogluon.tabular": "0.7.0",
    "scikit-learn": "1.2.2",
    "xgboost": "1.7.4",
}
SIMPLE_FEATURES = [
    "creatinine_min",
    "creatinine_max",
    "crea_divide_basecrea_min",
    "crea_divide_basecrea_max",
    "urineoutput_min",
    "urineoutput_mean",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument(
        "--allow-version-mismatch",
        action="store_true",
        help="Diagnostic override only; exact archived versions are required by default.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_environment(allow_mismatch: bool) -> dict[str, str]:
    observed = {
        "autogluon.tabular": metadata.version("autogluon.tabular"),
        "scikit-learn": sklearn.__version__,
        "xgboost": xgboost.__version__,
    }
    mismatches = {
        package: {"expected": expected, "observed": observed[package]}
        for package, expected in EXPECTED_VERSIONS.items()
        if observed[package] != expected
    }
    if mismatches and not allow_mismatch:
        raise RuntimeError(f"Archived environment mismatch: {mismatches}")
    return observed


def multiclass_brier(y: np.ndarray, prob: np.ndarray) -> float:
    onehot = np.eye(len(CLASSES))[y - 1]
    return float(np.mean(np.sum((prob - onehot) ** 2, axis=1)))


def calibration_parameters(y_binary: np.ndarray, prob: np.ndarray) -> tuple[float, float]:
    clipped = np.clip(prob, 1e-6, 1 - 1e-6)
    design = sm.add_constant(logit(clipped))
    try:
        fit = sm.GLM(y_binary, design, family=sm.families.Binomial()).fit()
        return float(fit.params[0]), float(fit.params[1])
    except Exception:
        return np.nan, np.nan


def expected_calibration_error(y_binary: np.ndarray, prob: np.ndarray) -> float:
    bins = pd.qcut(prob, q=10, labels=False, duplicates="drop")
    work = pd.DataFrame({"y": y_binary, "p": prob, "bin": bins})
    grouped = work.groupby("bin", observed=True)
    weights = grouped.size() / len(work)
    return float((weights * (grouped.y.mean() - grouped.p.mean()).abs()).sum())


def metric_values(y: np.ndarray, pred: np.ndarray, prob: np.ndarray) -> dict[str, float]:
    return {
        "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "precision_macro": float(
            precision_score(y, pred, average="macro", zero_division=0)
        ),
        "recall_macro": float(recall_score(y, pred, average="macro", zero_division=0)),
        "f1_macro": float(f1_score(y, pred, average="macro", zero_division=0)),
        "auc_ovo_macro": float(
            roc_auc_score(y, prob, multi_class="ovo", average="macro", labels=CLASSES)
        ),
        "auc_ovr_macro": float(
            roc_auc_score(y, prob, multi_class="ovr", average="macro", labels=CLASSES)
        ),
        "log_loss": float(log_loss(y, prob, labels=CLASSES)),
        "multiclass_brier": multiclass_brier(y, prob),
    }


def aggregate_evaluation(
    model_name: str,
    cohort: str,
    y: np.ndarray,
    pred: np.ndarray,
    prob: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    overall = pd.DataFrame(
        [{"model": model_name, "cohort": cohort, "n": len(y), **metric_values(y, pred, prob)}]
    )
    precision, recall, f1, support = precision_recall_fscore_support(
        y, pred, labels=CLASSES, zero_division=0
    )
    class_rows = []
    calibration_rows = []
    bin_rows = []
    for index, cls in enumerate(CLASSES):
        y_binary = (y == cls).astype(int)
        intercept, slope = calibration_parameters(y_binary, prob[:, index])
        class_rows.append(
            {
                "model": model_name,
                "cohort": cohort,
                "class": int(cls),
                "phenotype": PHENOTYPE[int(cls)],
                "support": int(support[index]),
                "prevalence": float(support[index] / len(y)),
                "precision": precision[index],
                "recall": recall[index],
                "f1": f1[index],
                "auc_ovr": roc_auc_score(y_binary, prob[:, index]),
            }
        )
        calibration_rows.append(
            {
                "model": model_name,
                "cohort": cohort,
                "class": int(cls),
                "phenotype": PHENOTYPE[int(cls)],
                "calibration_intercept": intercept,
                "calibration_slope": slope,
                "ece_deciles": expected_calibration_error(y_binary, prob[:, index]),
            }
        )
        bins = pd.qcut(prob[:, index], q=10, labels=False, duplicates="drop")
        work = pd.DataFrame({"y": y_binary, "p": prob[:, index], "bin": bins})
        grouped = work.groupby("bin", observed=True).agg(
            n=("y", "size"), mean_predicted=("p", "mean"), observed_frequency=("y", "mean")
        )
        for bin_index, row in grouped.reset_index().iterrows():
            bin_rows.append(
                {
                    "model": model_name,
                    "cohort": cohort,
                    "class": int(cls),
                    "phenotype": PHENOTYPE[int(cls)],
                    "bin": int(row["bin"]),
                    "n": int(row["n"]),
                    "mean_predicted": float(row["mean_predicted"]),
                    "observed_frequency": float(row["observed_frequency"]),
                }
            )

    pair_rows = []
    for a, b in [(1, 2), (1, 3), (2, 3)]:
        mask = np.isin(y, [a, b])
        score = prob[mask, a - 1] / (prob[mask, a - 1] + prob[mask, b - 1])
        pair_rows.append(
            {
                "model": model_name,
                "cohort": cohort,
                "comparison": f"{PHENOTYPE[a]} vs {PHENOTYPE[b]}",
                "n_pair": int(mask.sum()),
                "auc_pair_normalized": roc_auc_score((y[mask] == a).astype(int), score),
            }
        )
    return (
        overall,
        pd.DataFrame(class_rows),
        pd.DataFrame(calibration_rows),
        pd.DataFrame(bin_rows),
        pd.DataFrame(pair_rows),
    )


def paired_bootstrap_increment(
    cohort: str,
    y: np.ndarray,
    primary_pred: np.ndarray,
    primary_prob: np.ndarray,
    simple_pred: np.ndarray,
    simple_prob: np.ndarray,
    n_boot: int,
) -> pd.DataFrame:
    metric_names = ["balanced_accuracy", "f1_macro", "auc_ovo_macro", "multiclass_brier"]
    primary_point = metric_values(y, primary_pred, primary_prob)
    simple_point = metric_values(y, simple_pred, simple_prob)
    rng = np.random.default_rng(SEED + (1 if cohort.startswith("internal") else 2))
    by_class = {cls: np.flatnonzero(y == cls) for cls in CLASSES}
    draws = {metric: [] for metric in metric_names}
    for _ in range(n_boot):
        sampled = np.concatenate(
            [rng.choice(indices, size=len(indices), replace=True) for indices in by_class.values()]
        )
        rng.shuffle(sampled)
        primary = metric_values(y[sampled], primary_pred[sampled], primary_prob[sampled])
        simple = metric_values(y[sampled], simple_pred[sampled], simple_prob[sampled])
        for metric in metric_names:
            draws[metric].append(primary[metric] - simple[metric])
    rows = []
    for metric in metric_names:
        difference = primary_point[metric] - simple_point[metric]
        rows.append(
            {
                "cohort": cohort,
                "metric": metric,
                "archived_autogluon": primary_point[metric],
                "simple_logistic": simple_point[metric],
                "difference_autogluon_minus_simple": difference,
                "ci_low": np.quantile(draws[metric], 0.025),
                "ci_high": np.quantile(draws[metric], 0.975),
                "direction": (
                    "lower_is_better" if metric == "multiclass_brier" else "higher_is_better"
                ),
                "bootstrap_replicates": n_boot,
            }
        )
    return pd.DataFrame(rows)


def write_report(
    out_dir: Path,
    metrics: pd.DataFrame,
    calibration: pd.DataFrame,
    incremental: pd.DataFrame,
    versions: dict[str, str],
) -> None:
    def md(frame: pd.DataFrame) -> str:
        return frame.round(3).to_markdown(index=False)

    text = f"""# Archived AutoGluon classifier replay

## Status

The submitted `XGBoost_BAG_L2` model was successfully reloaded and evaluated using the
archived software versions: AutoGluon {versions['autogluon.tabular']}, scikit-learn
{versions['scikit-learn']}, and XGBoost {versions['xgboost']}. Only aggregate outputs
were exported; no patient identifiers or patient-level predictions were written.

## Primary model and simple comparator

{md(metrics)}

The comparator is a prespecified multinomial logistic regression using six transparent
first-24-hour kidney summaries: creatinine minimum/maximum, creatinine-to-baseline ratio
minimum/maximum, and urine-output minimum/mean. It is a deliberately simple reference,
not a proposed clinical score.

## One-versus-rest calibration

{md(calibration)}

Ideal calibration has intercept 0 and slope 1. ECE is an absolute decile-weighted
calibration error; lower is better.

## Paired incremental-value analysis

{md(incremental)}

Differences are AutoGluon minus simple logistic on the same test patients with
class-stratified bootstrap intervals. Positive differences favor AutoGluon for AUC,
F1, and balanced accuracy; negative differences favor AutoGluon for Brier score.
"""
    (out_dir / "ARCHIVED_AUTOGLUON_REPLAY.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    versions = validate_environment(args.allow_version_mismatch)

    from autogluon.tabular import TabularPredictor

    predictor = TabularPredictor.load(str(args.model_dir))
    input_dir = args.model_dir / "input"
    train = pd.read_csv(input_dir / "train_set.csv")
    tests = {
        "internal_MIMIC_eICU": pd.read_csv(input_dir / "test_set1.csv"),
        "external_AUMC": pd.read_csv(input_dir / "test_set2.csv"),
    }
    missing_simple = [column for column in SIMPLE_FEATURES if column not in train.columns]
    if missing_simple:
        raise ValueError(f"Missing simple-comparator features: {missing_simple}")

    simple_model = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    solver="lbfgs",
                    multi_class="multinomial",
                    C=1.0,
                    max_iter=5000,
                    random_state=SEED,
                ),
            ),
        ]
    )
    simple_model.fit(train[SIMPLE_FEATURES], train.groupHPD.astype(int))

    overall_parts = []
    class_parts = []
    calibration_parts = []
    calibration_bin_parts = []
    pairwise_parts = []
    incremental_parts = []
    for cohort, frame in tests.items():
        y = frame.groupHPD.astype(int).to_numpy()
        primary_df = predictor.predict_proba(frame, model="XGBoost_BAG_L2")
        primary_prob = primary_df.reindex(columns=[1.0, 2.0, 3.0]).to_numpy()
        primary_pred = CLASSES[np.argmax(primary_prob, axis=1)]
        simple_prob = simple_model.predict_proba(frame[SIMPLE_FEATURES])
        simple_pred = simple_model.predict(frame[SIMPLE_FEATURES]).astype(int)

        for model_name, pred, prob in [
            ("archived_AutoGluon_XGBoost_BAG_L2", primary_pred, primary_prob),
            ("simple_six_variable_logistic", simple_pred, simple_prob),
        ]:
            overall, per_class, calibration, calibration_bins, pairwise = aggregate_evaluation(
                model_name, cohort, y, pred, prob
            )
            overall_parts.append(overall)
            class_parts.append(per_class)
            calibration_parts.append(calibration)
            calibration_bin_parts.append(calibration_bins)
            pairwise_parts.append(pairwise)
        incremental_parts.append(
            paired_bootstrap_increment(
                cohort,
                y,
                primary_pred,
                primary_prob,
                simple_pred,
                simple_prob,
                args.bootstrap,
            )
        )

    overall_all = pd.concat(overall_parts, ignore_index=True)
    class_all = pd.concat(class_parts, ignore_index=True)
    calibration_all = pd.concat(calibration_parts, ignore_index=True)
    calibration_bins = pd.concat(calibration_bin_parts, ignore_index=True)
    pairwise = pd.concat(pairwise_parts, ignore_index=True)
    incremental = pd.concat(incremental_parts, ignore_index=True)

    overall_all.to_csv(args.out_dir / "primary_and_comparator_metrics.csv", index=False)
    class_all.to_csv(args.out_dir / "primary_and_comparator_class_metrics.csv", index=False)
    calibration_all.to_csv(args.out_dir / "primary_and_comparator_calibration.csv", index=False)
    calibration_bins.to_csv(args.out_dir / "primary_and_comparator_calibration_bins.csv", index=False)
    pairwise.to_csv(args.out_dir / "primary_and_comparator_pairwise_auc.csv", index=False)
    incremental.to_csv(args.out_dir / "paired_incremental_value.csv", index=False)
    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "versions": versions,
        "expected_versions": EXPECTED_VERSIONS,
        "environment_match": versions == EXPECTED_VERSIONS,
        "seed": SEED,
        "bootstrap_replicates": args.bootstrap,
        "simple_features": SIMPLE_FEATURES,
        "input_sha256": {
            filename: sha256(input_dir / filename)
            for filename in ["train_set.csv", "test_set1.csv", "test_set2.csv"]
        },
    }
    (args.out_dir / "archived_replay_environment.json").write_text(
        json.dumps(environment, indent=2), encoding="utf-8"
    )
    write_report(args.out_dir, overall_all, calibration_all, incremental, versions)
    print(overall_all.to_string(index=False))


if __name__ == "__main__":
    main()
