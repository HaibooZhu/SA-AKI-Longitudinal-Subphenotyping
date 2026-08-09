#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const workbookPaths = process.argv.slice(2);
if (workbookPaths.length === 0) {
  throw new Error("Pass at least one .xlsx path.");
}

const previewDir = process.env.JTIM_WORKBOOK_PREVIEW_DIR;
const jsonDir = process.env.JTIM_WORKBOOK_JSON_DIR;

for (const workbookPath of workbookPaths) {
  const blob = await FileBlob.load(workbookPath);
  const workbook = await SpreadsheetFile.importXlsx(blob);
  const workbookName = path.basename(workbookPath, path.extname(workbookPath));

  process.stdout.write(`\n## ${workbookName}\n`);
  const sheets = await workbook.inspect({
    kind: "sheet",
    include: "id,name",
    maxChars: 4000,
  });
  process.stdout.write(`${sheets.ndjson}\n`);

  for (const sheet of workbook.worksheets.items) {
    process.stdout.write(`\n### Sheet: ${sheet.name}\n`);
    const region = await workbook.inspect({
      kind: "region",
      sheetId: sheet.name,
      range: "A1:Z50",
      maxChars: 18000,
      tableMaxRows: 50,
      tableMaxCols: 26,
      tableMaxCellChars: 120,
    });
    process.stdout.write(`${region.ndjson}\n`);

    const formulas = await workbook.inspect({
      kind: "formula",
      sheetId: sheet.name,
      range: "A1:Z50",
      maxChars: 6000,
      options: { maxResults: 200 },
    });
    process.stdout.write(`FORMULAS\n${formulas.ndjson}\n`);

    if (jsonDir) {
      await fs.mkdir(jsonDir, { recursive: true });
      const safeSheetName = sheet.name.replace(/[^A-Za-z0-9_-]+/g, "_");
      const jsonPath = path.join(jsonDir, `${workbookName}_${safeSheetName}.json`);
      const values = sheet.getRange("B2:Z39").values;
      await fs.writeFile(
        jsonPath,
        `${JSON.stringify({ workbook: workbookName, sheet: sheet.name, range: "B2:Z39", values }, null, 2)}\n`,
      );
      process.stdout.write(`JSON ${jsonPath}\n`);
    }

    if (previewDir) {
      await fs.mkdir(previewDir, { recursive: true });
      const preview = await workbook.render({
        sheetName: sheet.name,
        autoCrop: "all",
        scale: 1,
        format: "png",
      });
      const safeSheetName = sheet.name.replace(/[^A-Za-z0-9_-]+/g, "_");
      const previewPath = path.join(previewDir, `${workbookName}_${safeSheetName}.png`);
      await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));
      process.stdout.write(`PREVIEW ${previewPath}\n`);
    }
  }
}
