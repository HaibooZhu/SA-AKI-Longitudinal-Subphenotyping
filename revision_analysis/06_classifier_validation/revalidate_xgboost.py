#!/usr/bin/env python3
"""Transparent classifier revalidation on the frozen manuscript splits.

This script does not claim bit-for-bit reproduction of the historical AutoGluon
stack. It (1) transcribes the archived AutoGluon XGBoost_BAG_L2 results and the
archived directed OvO ROC values, and (2) fits a fixed-parameter XGBoost model on
the same frozen train/test files as an independently rerunnable sensitivity
analysis.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import sklearn
import statsmodels.api as sm
import xgboost
from scipy.special import logit
from sklearn.calibration import calibration_curve
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    log_loss,
    precision_recall_fscore_support,
    precision_score,
    recall_score,
    roc_auc_score,
)
from xgboost import XGBClassifier


STYLE_DIR = Path(__file__).resolve().parents[1] / "07_tables_figures"
sys.path.insert(0, str(STYLE_DIR))
from publication_figure_style import (  # noqa: E402
    DOUBLE_COLUMN_IN,
    FONT_LEGEND,
    LINE_AUX,
    LINE_MAIN,
    PHENOTYPE_COLORS,
    add_panel_label,
    apply_publication_style,
    export_figure,
)


SEED = 1234
CLASSES = np.array([1, 2, 3], dtype=int)
PHENOTYPE = {1: "DR", 2: "RR", 3: "PW"}


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    frozen = (
        repo
        / "00_frozen_inputs/data_snapshot/remote_project_snapshot/07.autogluon/01.model"
        / "Result-a1234_selfv2_MimiceICU_AUMC_CorrMICfilt"
    )
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", type=Path, default=frozen)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=repo / "02_revision_outputs/reports/W4_classifier_validation",
    )
    parser.add_argument("--bootstrap", type=int, default=1000)
    return parser.parse_args()


def safe_auc(y: np.ndarray, score: np.ndarray) -> float:
    return float(roc_auc_score(y, score)) if np.unique(y).size == 2 else np.nan


def multiclass_brier(y: np.ndarray, prob: np.ndarray) -> float:
    onehot = np.eye(len(CLASSES))[y - 1]
    return float(np.mean(np.sum((prob - onehot) ** 2, axis=1)))


def expected_calibration_error(y_binary: np.ndarray, prob: np.ndarray) -> float:
    bins = pd.qcut(prob, q=10, labels=False, duplicates="drop")
    work = pd.DataFrame({"y": y_binary, "p": prob, "bin": bins})
    grouped = work.groupby("bin", observed=True)
    weights = grouped.size() / len(work)
    return float((weights * (grouped.y.mean() - grouped.p.mean()).abs()).sum())


def calibration_parameters(y_binary: np.ndarray, prob: np.ndarray) -> tuple[float, float]:
    clipped = np.clip(prob, 1e-6, 1 - 1e-6)
    design = sm.add_constant(logit(clipped))
    try:
        fit = sm.GLM(y_binary, design, family=sm.families.Binomial()).fit()
        return float(fit.params[0]), float(fit.params[1])
    except Exception:
        return np.nan, np.nan


def performance_tables(
    cohort: str, y: np.ndarray, pred: np.ndarray, prob: np.ndarray
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    overall = pd.DataFrame(
        [
            {
                "cohort": cohort,
                "n": len(y),
                "accuracy": accuracy_score(y, pred),
                "balanced_accuracy": balanced_accuracy_score(y, pred),
                "precision_macro": precision_score(y, pred, average="macro", zero_division=0),
                "recall_macro": recall_score(y, pred, average="macro", zero_division=0),
                "f1_macro": f1_score(y, pred, average="macro", zero_division=0),
                "auc_ovo_macro": roc_auc_score(y, prob, multi_class="ovo", average="macro", labels=CLASSES),
                "auc_ovr_macro": roc_auc_score(y, prob, multi_class="ovr", average="macro", labels=CLASSES),
                "log_loss": log_loss(y, prob, labels=CLASSES),
                "multiclass_brier": multiclass_brier(y, prob),
            }
        ]
    )

    cm = confusion_matrix(y, pred, labels=CLASSES)
    precision, recall, f1, support = precision_recall_fscore_support(
        y, pred, labels=CLASSES, zero_division=0
    )
    class_rows = []
    calibration_rows = []
    for idx, cls in enumerate(CLASSES):
        tp = cm[idx, idx]
        fn = cm[idx, :].sum() - tp
        fp = cm[:, idx].sum() - tp
        tn = cm.sum() - tp - fn - fp
        y_binary = (y == cls).astype(int)
        intercept, slope = calibration_parameters(y_binary, prob[:, idx])
        class_rows.append(
            {
                "cohort": cohort,
                "class": int(cls),
                "phenotype": PHENOTYPE[int(cls)],
                "support": int(support[idx]),
                "prevalence": float(support[idx] / len(y)),
                "precision": precision[idx],
                "sensitivity_recall": recall[idx],
                "specificity": tn / (tn + fp),
                "npv": tn / (tn + fn),
                "f1": f1[idx],
                "auc_ovr": safe_auc(y_binary, prob[:, idx]),
            }
        )
        calibration_rows.append(
            {
                "cohort": cohort,
                "class": int(cls),
                "phenotype": PHENOTYPE[int(cls)],
                "calibration_intercept": intercept,
                "calibration_slope": slope,
                "ece_deciles": expected_calibration_error(y_binary, prob[:, idx]),
            }
        )

    pair_rows = []
    for positive in CLASSES:
        for comparator in CLASSES:
            if positive == comparator:
                continue
            mask = np.isin(y, [positive, comparator])
            y_binary = (y[mask] == positive).astype(int)
            positive_idx = int(np.where(CLASSES == positive)[0][0])
            comparator_idx = int(np.where(CLASSES == comparator)[0][0])
            raw_auc = safe_auc(y_binary, prob[mask, positive_idx])
            denom = prob[mask, positive_idx] + prob[mask, comparator_idx]
            normalized = np.divide(
                prob[mask, positive_idx],
                denom,
                out=np.full(mask.sum(), 0.5),
                where=denom > 0,
            )
            pair_rows.append(
                {
                    "cohort": cohort,
                    "positive_class": int(positive),
                    "positive_phenotype": PHENOTYPE[int(positive)],
                    "comparator_class": int(comparator),
                    "comparator_phenotype": PHENOTYPE[int(comparator)],
                    "n_pair": int(mask.sum()),
                    "auc_directed_raw_probability": raw_auc,
                    "auc_pair_normalized": safe_auc(y_binary, normalized),
                }
            )
    return overall, pd.DataFrame(class_rows), pd.DataFrame(pair_rows), pd.DataFrame(calibration_rows)


def bootstrap_key_metrics(
    cohort: str, y: np.ndarray, pred: np.ndarray, prob: np.ndarray, n_boot: int
) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    by_class = {cls: np.flatnonzero(y == cls) for cls in CLASSES}
    rows = []
    for _ in range(n_boot):
        sampled = np.concatenate(
            [rng.choice(idx, size=len(idx), replace=True) for idx in by_class.values()]
        )
        rng.shuffle(sampled)
        yb, predb, probb = y[sampled], pred[sampled], prob[sampled]
        row = {
            "accuracy": accuracy_score(yb, predb),
            "balanced_accuracy": balanced_accuracy_score(yb, predb),
            "f1_macro": f1_score(yb, predb, average="macro", zero_division=0),
            "auc_ovo_macro": roc_auc_score(
                yb, probb, multi_class="ovo", average="macro", labels=CLASSES
            ),
        }
        for a, b in [(1, 2), (1, 3), (2, 3)]:
            mask = np.isin(yb, [a, b])
            ia, ib = a - 1, b - 1
            score = probb[mask, ia] / (probb[mask, ia] + probb[mask, ib])
            row[f"auc_{PHENOTYPE[a]}_vs_{PHENOTYPE[b]}"] = safe_auc(
                (yb[mask] == a).astype(int), score
            )
        rows.append(row)
    boot = pd.DataFrame(rows)
    return pd.DataFrame(
        [
            {
                "cohort": cohort,
                "metric": metric,
                "estimate": boot[metric].median(),
                "ci_low": boot[metric].quantile(0.025),
                "ci_high": boot[metric].quantile(0.975),
                "bootstrap_replicates": n_boot,
            }
            for metric in boot.columns
        ]
    )


def original_tables(model_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    result_dir = model_dir / "result"
    sources = {
        "internal_MIMIC_eICU": result_dir / "result_test_['eicu', 'mimic'].csv",
        "external_AUMC": result_dir / "result_test_aumcdb.csv",
    }
    rows = []
    for cohort, path in sources.items():
        source = pd.read_csv(path)
        row = source.loc[source.model == "XGBoost_BAG_L2"].copy()
        if len(row) != 1:
            raise ValueError(f"Expected one XGBoost_BAG_L2 row in {path}")
        row.insert(0, "cohort", cohort)
        row.insert(1, "source_file", path.name)
        rows.append(row)
    performance = pd.concat(rows, ignore_index=True).drop(columns=["Unnamed: 0"])

    archived_auc = {
        "internal_MIMIC_eICU": {
            ("DR", "RR"): 0.8505,
            ("DR", "PW"): 0.6168,
            ("RR", "DR"): 0.8635,
            ("RR", "PW"): 0.9427,
            ("PW", "DR"): 0.7998,
            ("PW", "RR"): 0.9426,
        },
        "external_AUMC": {
            ("DR", "RR"): 0.6680,
            ("DR", "PW"): 0.5931,
            ("RR", "DR"): 0.8371,
            ("RR", "PW"): 0.9245,
            ("PW", "DR"): 0.6730,
            ("PW", "RR"): 0.9187,
        },
    }
    pairwise = pd.DataFrame(
        [
            {
                "cohort": cohort,
                "positive_phenotype": pair[0],
                "comparator_phenotype": pair[1],
                "auc_directed_raw_probability": auc,
                "source_file": (
                    "all_OVO_ROC_Curve_on_MIMIC_eICU_test_set.pdf"
                    if cohort == "internal_MIMIC_eICU"
                    else "all_OVO_ROC_Curve_on_AUMCdb_test_set.pdf"
                ),
            }
            for cohort, values in archived_auc.items()
            for pair, auc in values.items()
        ]
    )
    return performance, pairwise


def plot_validation(
    evaluations: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], out_dir: Path
) -> None:
    apply_publication_style()
    colors = {group: PHENOTYPE_COLORS[name] for group, name in PHENOTYPE.items()}
    cohort_titles = {
        "internal_MIMIC_eICU": "Internal: MIMIC-IV + eICU-CRD",
        "external_AUMC": "External: AmsterdamUMCdb",
    }
    pale_blues = LinearSegmentedColormap.from_list(
        "jtim_classifier_blues", ["#F7F9FC", "#D5E4F0", "#8DB9D8"]
    )
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(DOUBLE_COLUMN_IN, 5.4),
        constrained_layout=True,
        gridspec_kw={"height_ratios": [0.82, 1.18]},
    )
    for col, (cohort, (y, pred, prob)) in enumerate(evaluations.items()):
        cm = confusion_matrix(y, pred, labels=CLASSES, normalize="true")
        ax = axes[0, col]
        image = ax.imshow(cm, vmin=0, vmax=1, cmap=pale_blues)
        for row in range(3):
            for column in range(3):
                color = "white" if cm[row, column] > 0.78 else "#303030"
                ax.text(column, row, f"{cm[row, column]:.2f}", ha="center", va="center", color=color)
        ax.set_xticks(range(3), [PHENOTYPE[x] for x in CLASSES])
        ax.set_yticks(range(3), [PHENOTYPE[x] for x in CLASSES])
        ax.set_xlabel("Predicted phenotype")
        ax.set_ylabel("Observed phenotype")
        ax.set_title(cohort_titles[cohort], fontweight="bold", pad=5)
        ax.spines[:].set_visible(False)

        ax = axes[1, col]
        ax.plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=LINE_AUX, label="Ideal")
        for idx, cls in enumerate(CLASSES):
            observed, predicted = calibration_curve(
                (y == cls).astype(int), prob[:, idx], n_bins=10, strategy="quantile"
            )
            ax.plot(
                predicted,
                observed,
                marker="o",
                markersize=3.0,
                linewidth=LINE_MAIN,
                color=colors[int(cls)],
                label=PHENOTYPE[int(cls)],
            )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Predicted probability")
        ax.set_ylabel("Observed frequency")
        ax.set_title("Calibration", pad=5)
    colorbar = fig.colorbar(
        image,
        ax=axes[0, :],
        shrink=0.68,
        fraction=0.032,
        pad=0.025,
        label="Row-normalized proportion",
    )
    colorbar.outline.set_linewidth(0.6)
    for label, ax in zip(["a", "b", "c", "d"], axes.ravel()):
        add_panel_label(ax, label)
    shared_handles = [
        Line2D([], [], color="#777777", linestyle="--", lw=LINE_AUX, label="Ideal"),
        *[
            Line2D([], [], color=colors[cls], marker="o", markersize=3, lw=LINE_MAIN, label=PHENOTYPE[cls])
            for cls in CLASSES
        ],
    ]
    fig.legend(
        handles=shared_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.505),
        ncol=4,
        fontsize=FONT_LEGEND,
        columnspacing=1.0,
        handletextpad=0.35,
    )
    stem = out_dir / "W4_classifier_revalidation"
    export_figure(fig, stem)
    plt.close(fig)


def plot_archived_calibration_bins(bins: pd.DataFrame, out_dir: Path) -> None:
    """Plot aggregate calibration bins exported by the exact-version replay."""
    apply_publication_style()
    colors = PHENOTYPE_COLORS
    models = [
        ("archived_AutoGluon_XGBoost_BAG_L2", "Archived AutoGluon XGBoost"),
        ("simple_six_variable_logistic", "Six-variable logistic comparator"),
    ]
    cohorts = [
        ("internal_MIMIC_eICU", "Internal MIMIC-IV/eICU"),
        ("external_AUMC", "External AUMC"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COLUMN_IN, 5.8), constrained_layout=True)
    for row, (model, model_label) in enumerate(models):
        for column, (cohort, cohort_label) in enumerate(cohorts):
            ax = axes[row, column]
            ax.plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=LINE_AUX)
            subset = bins.loc[bins.model.eq(model) & bins.cohort.eq(cohort)]
            for phenotype in ["DR", "RR", "PW"]:
                values = subset.loc[subset.phenotype.eq(phenotype)].sort_values("bin")
                ax.plot(
                    values.mean_predicted,
                    values.observed_frequency,
                    marker="o",
                    markersize=3.0,
                    linewidth=LINE_MAIN,
                    color=colors[phenotype],
                    label=phenotype,
                )
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_xlabel("Predicted probability")
            ax.set_ylabel("Observed frequency")
            ax.set_title(f"{model_label}\n{cohort_label}")
    for label, ax in zip(["a", "b", "c", "d"], axes.ravel()):
        add_panel_label(ax, label)
    handles = [
        Line2D([], [], color="#777777", linestyle="--", lw=LINE_AUX, label="Ideal"),
        *[
            Line2D([], [], color=colors[name], marker="o", markersize=3, lw=LINE_MAIN, label=name)
            for name in ["DR", "RR", "PW"]
        ],
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.015), ncol=4)
    stem = out_dir / "W4_archived_model_calibration_comparator"
    export_figure(fig, stem)
    plt.close(fig)


def write_report(
    out_dir: Path,
    original_performance: pd.DataFrame,
    original_pairwise: pd.DataFrame,
    overall: pd.DataFrame,
    per_class: pd.DataFrame,
    bootstrap: pd.DataFrame,
    replay_metrics: pd.DataFrame | None = None,
    replay_calibration: pd.DataFrame | None = None,
    incremental: pd.DataFrame | None = None,
) -> None:
    def md(df: pd.DataFrame, digits: int = 3) -> str:
        return df.round(digits).to_markdown(index=False)

    original_cols = [
        "cohort",
        "OVO_macro_AUC_test",
        "accuracy",
        "precision_macro",
        "recall_macro",
        "f1_macro",
        "precision_C1",
        "precision_C2",
        "precision_C3",
        "recall_C1",
        "recall_C2",
        "recall_C3",
        "f1_C1",
        "f1_C2",
        "f1_C3",
    ]
    if replay_metrics is not None:
        replay_section = f"""
## Exact-version replay of the submitted AutoGluon model

The archived `XGBoost_BAG_L2` predictor was reloaded in AutoGluon 0.7.0,
scikit-learn 1.2.2, and XGBoost 1.7.4. The comparator is a fixed multinomial
logistic regression using six transparent first-24-hour kidney summaries:
creatinine minimum/maximum, creatinine-to-baseline ratio minimum/maximum, and
urine-output minimum/mean.

{md(replay_metrics)}

### Calibration of the submitted model and comparator

{md(replay_calibration)}

Ideal calibration has intercept 0 and slope 1; lower ECE is better. The submitted
model is reasonably calibrated on the internal test set but markedly miscalibrated
in AUMC, especially for DR and RR.

### Paired incremental value

{md(incremental)}

Differences are archived AutoGluon minus simple logistic on the same patients with
class-stratified bootstrap intervals. The archived model improves all four metrics
internally. In external AUMC, it does not improve macro OvO AUC and performs worse on
balanced accuracy, macro F1, and Brier score. Consequently, incremental clinical value
beyond the simple kidney summary model is not established externally.
"""
    else:
        replay_section = """
## Exact-version replay status

Exact-version AutoGluon replay outputs were not available in this run. Do not make a
calibration or incremental-value claim from the fixed-parameter sensitivity model.
"""
    text = f"""# W4 Classifier validation audit

## Bottom line

The archived results support moderate internal discrimination with degradation in the
external AUMC cohort. They do not support presenting a single AUC of 0.943 as the
overall classifier performance. In particular, the archived directed DR-versus-PW
AUC was 0.617 internally and 0.593 externally, while the macro OvO AUC decreased
from 0.836 to 0.769. The classifier should therefore be described as a retrospective
research tool requiring prospective validation, not as an actionable bedside tool.

## Archived AutoGluon XGBoost_BAG_L2 performance

{md(original_performance[original_cols])}

### Complete archived directed OvO AUCs

{md(original_pairwise.drop(columns='source_file'))}

These six directed values use each named positive class's raw multiclass probability;
the two directions of a class pair are therefore not mathematically constrained to be
identical. Values were transcribed from the archived vector PDFs and visually checked
against rendered pages.

{replay_section}

## Independent fixed-parameter XGBoost revalidation

This is a transparent sensitivity analysis on the same frozen splits, not a bit-for-bit
reproduction of the historical AutoGluon two-layer bagged stack. Patient identifiers
were excluded; all preprocessing was frozen before model fitting; the training file was
used for fitting and the two test files were not used for tuning.

{md(overall)}

### Per-class performance

{md(per_class)}

### Stratified bootstrap intervals

{md(bootstrap)}

## Reporting decision

1. Foreground macro OvO AUC and macro F1 for internal and external validation.
2. Report all directed pairwise OvO AUCs and all class-specific precision, recall, and
   F1 values rather than selecting the most favorable pair.
3. Add confusion matrices and calibration plots to the supplement.
4. Remove claims of immediate clinical deployment, treatment assignment, or reliable
   bedside differentiation between all three phenotypes.
5. State that external incremental value over the simple comparator was not shown and
   that prospective validation and recalibration are required.
"""
    (out_dir / "W4_CLASSIFIER_VALIDATION.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    input_dir = args.model_dir / "input"

    train = pd.read_csv(input_dir / "train_set.csv")
    tests = {
        "internal_MIMIC_eICU": pd.read_csv(input_dir / "test_set1.csv"),
        "external_AUMC": pd.read_csv(input_dir / "test_set2.csv"),
    }
    feature_cols = [c for c in train.columns if c not in {"stay_id", "groupHPD"}]
    y_train = train.groupHPD.astype(int).to_numpy()

    model = XGBClassifier(
        objective="multi:softprob",
        num_class=3,
        n_estimators=300,
        learning_rate=0.1,
        max_depth=6,
        min_child_weight=1,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        eval_metric="mlogloss",
        random_state=SEED,
        n_jobs=-1,
        tree_method="hist",
    )
    model.fit(train[feature_cols], y_train - 1)
    model.save_model(args.out_dir / "independent_xgboost_model.json")

    evaluations: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    overall_parts, class_parts, pair_parts, calibration_parts, bootstrap_parts = [], [], [], [], []
    for cohort, frame in tests.items():
        y = frame.groupHPD.astype(int).to_numpy()
        prob = model.predict_proba(frame[feature_cols])
        pred = CLASSES[np.argmax(prob, axis=1)]
        evaluations[cohort] = (y, pred, prob)
        overall, per_class, pairwise, calibration = performance_tables(cohort, y, pred, prob)
        overall_parts.append(overall)
        class_parts.append(per_class)
        pair_parts.append(pairwise)
        calibration_parts.append(calibration)
        bootstrap_parts.append(bootstrap_key_metrics(cohort, y, pred, prob, args.bootstrap))

    original_performance, original_pairwise = original_tables(args.model_dir)
    overall_all = pd.concat(overall_parts, ignore_index=True)
    class_all = pd.concat(class_parts, ignore_index=True)
    pair_all = pd.concat(pair_parts, ignore_index=True)
    calibration_all = pd.concat(calibration_parts, ignore_index=True)
    bootstrap_all = pd.concat(bootstrap_parts, ignore_index=True)

    original_performance.to_csv(args.out_dir / "archived_xgboost_bag_l2_metrics.csv", index=False)
    original_pairwise.to_csv(args.out_dir / "archived_directed_pairwise_auc.csv", index=False)
    overall_all.to_csv(args.out_dir / "revalidation_overall_metrics.csv", index=False)
    class_all.to_csv(args.out_dir / "revalidation_class_metrics.csv", index=False)
    pair_all.to_csv(args.out_dir / "revalidation_pairwise_auc.csv", index=False)
    calibration_all.to_csv(args.out_dir / "revalidation_calibration.csv", index=False)
    bootstrap_all.to_csv(args.out_dir / "revalidation_bootstrap_ci.csv", index=False)

    plot_validation(evaluations, args.out_dir)
    replay_metrics_path = args.out_dir / "primary_and_comparator_metrics.csv"
    replay_calibration_path = args.out_dir / "primary_and_comparator_calibration.csv"
    replay_bins_path = args.out_dir / "primary_and_comparator_calibration_bins.csv"
    incremental_path = args.out_dir / "paired_incremental_value.csv"
    if all(
        path.exists()
        for path in [
            replay_metrics_path,
            replay_calibration_path,
            replay_bins_path,
            incremental_path,
        ]
    ):
        replay_metrics = pd.read_csv(replay_metrics_path)
        replay_calibration = pd.read_csv(replay_calibration_path)
        replay_bins = pd.read_csv(replay_bins_path)
        incremental = pd.read_csv(incremental_path)
        plot_archived_calibration_bins(replay_bins, args.out_dir)
    else:
        replay_metrics = replay_calibration = incremental = None

    write_report(
        args.out_dir,
        original_performance,
        original_pairwise,
        overall_all,
        class_all,
        bootstrap_all,
        replay_metrics=replay_metrics,
        replay_calibration=replay_calibration,
        incremental=incremental,
    )
    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
        "statsmodels": sm.__version__,
        "xgboost": xgboost.__version__,
        "seed": SEED,
        "n_bootstrap": args.bootstrap,
        "n_features": len(feature_cols),
        "features": feature_cols,
    }
    (args.out_dir / "environment_and_features.json").write_text(
        json.dumps(environment, indent=2), encoding="utf-8"
    )
    print(overall_all.to_string(index=False))


if __name__ == "__main__":
    main()
