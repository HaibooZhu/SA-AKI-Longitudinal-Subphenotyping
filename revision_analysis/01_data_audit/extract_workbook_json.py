#!/usr/bin/env python3
"""Extract audited workbook cells to the JSON format used by the Table S3 audit."""

from __future__ import annotations

import argparse
import json
import re
from datetime import date, datetime, time
from pathlib import Path
from typing import Any

from openpyxl import load_workbook


AUDIT_RANGE = "B2:Z39"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract B2:Z39 from every worksheet using cached workbook values."
    )
    parser.add_argument("workbooks", nargs="+", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def json_value(value: Any) -> Any:
    if isinstance(value, (datetime, date, time)):
        return value.isoformat()
    return value


def safe_sheet_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]+", "_", name)


def extract_workbook(path: Path, output_dir: Path) -> None:
    workbook = load_workbook(path, read_only=True, data_only=True)
    workbook_name = path.stem
    try:
        for worksheet in workbook.worksheets:
            values = [
                [json_value(cell.value) for cell in row]
                for row in worksheet[AUDIT_RANGE]
            ]
            payload = {
                "workbook": workbook_name,
                "sheet": worksheet.title,
                "range": AUDIT_RANGE,
                "values": values,
            }
            output_path = output_dir / (
                f"{workbook_name}_{safe_sheet_name(worksheet.title)}.json"
            )
            output_path.write_text(
                json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
            print(output_path)
    finally:
        workbook.close()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for workbook_path in args.workbooks:
        extract_workbook(workbook_path.resolve(strict=True), output_dir)


if __name__ == "__main__":
    main()
