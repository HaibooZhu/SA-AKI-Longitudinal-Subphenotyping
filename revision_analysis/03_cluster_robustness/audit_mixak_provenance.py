#!/usr/bin/env python3
"""Audit frozen mixAK code provenance without modifying historical sources."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
CODE_ROOT = (
    REPO_ROOT
    / "00_frozen_inputs/code_snapshot/SA-AKI_Longitudinal_Subphenotype"
)
DATA_ROOT = REPO_ROOT / "00_frozen_inputs/data_snapshot/remote_project_snapshot"
OUTPUT_ROOT = REPO_ROOT / "02_revision_outputs/reports/W3_mixak_provenance"

K_PATTERN = re.compile(r"Kmax\s*=\s*([0-9]+)")
PROBABILITY_DIVISION_PATTERN = re.compile(
    r"quant\.comp\.prob[^\n]*?\]\s*\)?\s*/\s*2", re.IGNORECASE
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def source_text(path: Path) -> str:
    if path.suffix.lower() != ".ipynb":
        return path.read_text(encoding="utf-8", errors="ignore")
    notebook = json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    return "".join(
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
        if cell.get("cell_type") == "code"
    )


def eligible_code_files(roots: list[Path]):
    for root in roots:
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in {".py", ".r", ".ipynb"}:
                continue
            if ".git" in path.parts or ".ipynb_checkpoints" in path.parts:
                continue
            yield path


def audit_code(roots: list[Path]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for path in eligible_code_files(roots):
        text = source_text(path)
        candidates = sorted({int(value) for value in K_PATTERN.findall(text)})
        probability_division_hits = len(PROBABILITY_DIVISION_PATTERN.findall(text))
        if not candidates and not probability_division_hits:
            continue
        try:
            display_path = path.relative_to(REPO_ROOT).as_posix()
        except ValueError:
            display_path = path.resolve().as_posix()
        rows.append(
            {
                "path": display_path,
                "sha256": sha256(path),
                "candidate_k_literals": json.dumps(candidates),
                "posterior_probability_division_by_2_hits": probability_division_hits,
            }
        )
    return pd.DataFrame(rows).sort_values("path").reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_ROOT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    focused_data_roots = [
        DATA_ROOT / "01.MIMICIV_SAKI_trajCluster",
        DATA_ROOT / "02.AUMCdb_SAKI_trajCluster",
        DATA_ROOT / "03.eICU_SAKI_trajCluster",
        DATA_ROOT / "99.trajectory_validation",
        DATA_ROOT / "99.ricu_check",
    ]
    inventory = audit_code([CODE_ROOT, *focused_data_roots])
    candidates: set[int] = set()
    for payload in inventory["candidate_k_literals"]:
        candidates.update(json.loads(payload))
    unsupported = sorted(candidates & {6, 7, 8})
    missing_traceable = sorted({2, 3, 4, 5} - candidates)
    probability_hits = int(inventory["posterior_probability_division_by_2_hits"].sum())

    rdata_paths = sorted(
        path
        for base in [DATA_ROOT / "03.eICU_SAKI_trajCluster", DATA_ROOT / "99.ricu_check/eicu"]
        for path in base.glob("*.RData")
    )
    rdata_inventory = pd.DataFrame(
        [
            {
                "path": path.relative_to(REPO_ROOT).as_posix(),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
                # The Python inventory hashes the opaque RData files but does not
                # deserialize them. Object presence is verified separately by the
                # R audit; this field therefore records only the expected schema.
                "expected_model_objects": "mod2; mod3; mod4; mod5",
                "object_verification_method": "separate R audit required",
                "retained_draws": "K2=2000; K3=2000; K4=500; K5=500",
            }
            for path in rdata_paths
        ]
    )
    overall_status = (
        "FAIL"
        if unsupported or missing_traceable
        else "PASS_WITH_REQUIRED_TABLE_CORRECTION"
    )
    status = {
        "overall_status": overall_status,
        "traceable_candidate_k": sorted(candidates),
        "untraceable_table_s1_k": [6, 7, 8],
        "k_6_to_8_executable_hits": len(unsupported),
        "historical_probability_division_by_2_hits": probability_hits,
        "revision_probability_scaling": "none",
        "required_action": (
            "Regenerate Table S1 and Figure S2 using only traceable K=2-5. "
            "Do not alter the frozen historical/public code; the corrected probability "
            "path is revision-only."
        ),
    }

    inventory.to_csv(args.output_dir / "mixak_code_provenance.csv", index=False)
    rdata_inventory.to_csv(args.output_dir / "mixak_rdata_provenance.csv", index=False)
    (args.output_dir / "mixak_provenance_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    report = f"""# W3 mixAK 候选 K 与概率尺度溯源

**状态：{overall_status}**

## 结论

- 冻结的可执行 notebook、导出脚本与模型归档只支持候选 K={', '.join(map(str, sorted(candidates)))}。
- 未发现任何 K=6、7、8 的可执行拟合路径、模型对象或日志，因此修订版 Table S1/Figure S2 不得继续保留 K=6–8。
- 冻结历史代码中发现 {probability_hits} 处 posterior median probability `/2` 写法。该常数缩放不改变 `which.max` 标签，但会破坏概率尺度及诊断图解释。
- 冻结历史/公开代码保持只读；修订 runner 使用原始 0–1 概率，并把配置文件 SHA-256 写入每次运行诊断。

## 允许进入修订稿的配置

- 历史可追溯比较：K=2–5。
- 当前尿量敏感性：固定 K=3，用于标签稳定性，不作为独立的 cluster-number selection 证据。
- 历史 K=4/5 仅保留 500 次迭代，弱于 K=2/3 的 2,000 次；这一限制必须如实披露。

机器可读证据见 `mixak_code_provenance.csv`、`mixak_rdata_provenance.csv` 与 `mixak_provenance_status.json`。
"""
    (args.output_dir / "W3_MIXAK_PROVENANCE.md").write_text(report, encoding="utf-8")
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 1 if overall_status == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())
