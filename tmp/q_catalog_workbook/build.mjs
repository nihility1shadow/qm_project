import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = "E:/project_qm";
const catalogPath = `${root}/scan_cloud_8m_t150_20260809/all_parameter_q_catalog.csv`;
const outputPath = `${root}/scan_cloud_8m_t150_20260809/parameter_q_lookup.xlsx`;
const previewDir = `${root}/scan_cloud_8m_t150_20260809/workbook_previews`;

const csvText = await fs.readFile(catalogPath, "utf8");
const workbook = await Workbook.fromCSV(csvText, { sheetName: "All Evaluations" });
const all = workbook.worksheets.getItem("All Evaluations");
const used = all.getUsedRange();
const values = used.values;
const headers = values[0];
const rows = values.slice(1);
const lastRow = rows.length + 1;
const lastColumn = headers.length;
const uniqueParameterCount = new Set(rows.map((row) => JSON.stringify([
  row[1], row[7], row[8], row[10], row[12], row[13],
]))).size;

const currentRows = rows
  .filter((row) => Number(row[10]) === 8000000
    && Number(row[12]) === 150
    && Number(row[13]) === 130)
  .sort((left, right) => Number(right[14]) - Number(left[14]));
const validatedRows = currentRows
  .filter((row) => Number(row[11]) >= 6)
  .sort((left, right) => Number(right[14]) - Number(left[14]));

const current = workbook.worksheets.add("8M 150t");
current.getRangeByIndexes(0, 0, currentRows.length + 1, lastColumn).values = [
  headers,
  ...currentRows,
];

const summary = workbook.worksheets.add("Summary");
summary.getRange("A1:H1").merge();
summary.getRange("A1").values = [["Poisson parameter convergence catalog"]];
summary.getRange("A3:B11").values = [
  ["Metric", "Value"],
  ["All evaluations", null],
  ["Unique parameter settings", uniqueParameterCount],
  ["8M / 150 a.u. evaluations", null],
  ["Target passes (Q > 4 through 130)", null],
  ["Best exploratory strict Q", null],
  ["Best validated strict Q", validatedRows.length ? Number(validatedRows[0][14]) : null],
  ["Best validated case", validatedRows.length ? validatedRows[0][6] : "none"],
  ["Catalog updated", new Date()],
];
summary.getRange("B4").formulas = [[`=COUNTA('All Evaluations'!$A$2:$A$${lastRow})`]];
const currentLastRow = Math.max(currentRows.length + 1, 2);
summary.getRange("B6").formulas = [[`=COUNTA('8M 150t'!$A$2:$A$${currentLastRow})`]];
summary.getRange("B7").formulas = [[`=COUNTIF('8M 150t'!$Q$2:$Q$${currentLastRow},"TRUE")`]];
summary.getRange("B8").formulas = [[`=MAX('8M 150t'!$O$2:$O$${currentLastRow})`]];
summary.getRange("B11").format.numberFormat = "yyyy-mm-dd hh:mm";

summary.getRange("D3:H10").values = [
  ["Search rule", "Setting", null, null, null],
  ["Paths per run", 8000000, null, null, null],
  ["Exploratory repeats", 3, null, null, null],
  ["Final validation repeats", 6, null, null, null],
  ["Simulation end", 150, "a.u.", null, null],
  ["Required interval", 130, "a.u.", null, null],
  ["Pass threshold", 4, "minimum active-orbital Q", null, null],
  ["Round retention", 0.3, "top fraction", null, null],
];
summary.getRange("E10").format.numberFormat = "0%";

const headerStyle = {
  fill: "#0F4C5C",
  font: { bold: true, color: "#FFFFFF" },
  verticalAlignment: "center",
};
const titleStyle = {
  fill: "#17324D",
  font: { bold: true, color: "#FFFFFF", fontSize: 18 },
  horizontalAlignment: "center",
  verticalAlignment: "center",
};

summary.showGridLines = false;
summary.getRange("A1:H1").format = titleStyle;
summary.getRange("A1:H1").format.rowHeight = 34;
summary.getRange("A3:B3").format = headerStyle;
summary.getRange("D3:H3").format = headerStyle;
summary.getRange("A3:B11").format.borders = {
  insideHorizontal: { style: "thin", color: "#CBD5E1" },
};
summary.getRange("D3:H10").format.borders = {
  insideHorizontal: { style: "thin", color: "#CBD5E1" },
};
summary.getRange("A1:A10").format.columnWidth = 34;
summary.getRange("B1:B10").format.columnWidth = 23;
summary.getRange("C1:C10").format.columnWidth = 4;
summary.getRange("D1:D10").format.columnWidth = 24;
summary.getRange("E1:E10").format.columnWidth = 16;
summary.getRange("F1:H10").format.columnWidth = 19;
summary.freezePanes.freezeRows(1);

function styleDataSheet(sheet, rowCount, tableName) {
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(7);
  const header = sheet.getRangeByIndexes(0, 0, 1, lastColumn);
  header.format = headerStyle;
  header.format.rowHeight = 30;
  header.format.wrapText = true;
  sheet.getRange(`A1:A${rowCount}`).format.columnWidth = 28;
  sheet.getRange(`B1:B${rowCount}`).format.columnWidth = 18;
  sheet.getRange(`C1:C${rowCount}`).format.columnWidth = 48;
  sheet.getRange(`D1:E${rowCount}`).format.columnWidth = 12;
  sheet.getRange(`F1:F${rowCount}`).format.columnWidth = 16;
  sheet.getRange(`G1:G${rowCount}`).format.columnWidth = 30;
  sheet.getRange(`H1:I${rowCount}`).format.columnWidth = 12;
  sheet.getRange(`J1:J${rowCount}`).format.columnWidth = 15;
  sheet.getRange(`K1:Q${rowCount}`).format.columnWidth = 14;
  sheet.getRange(`R1:S${rowCount}`).format.columnWidth = 18;
  sheet.getRangeByIndexes(0, 19, rowCount, lastColumn - 19).format.columnWidth = 15;
  if (rowCount > 1) {
    sheet.getRange(`H2:I${rowCount}`).format.numberFormat = "0.0000E+00";
    sheet.getRange(`K2:L${rowCount}`).format.numberFormat = "#,##0";
    sheet.getRange(`M2:N${rowCount}`).format.numberFormat = "0.0";
    sheet.getRange(`O2:O${rowCount}`).format.numberFormat = "0.000";
    sheet.getRange(`T2:AN${rowCount}`).format.numberFormat = "0.000E+00";
    sheet.getRange(`O2:O${rowCount}`).conditionalFormats.add("colorScale", {
      thresholds: ["min", "50%", "max"],
      colors: ["#FEE2E2", "#FEF3C7", "#DCFCE7"],
    });
    sheet.getRange(`F2:F${rowCount}`).conditionalFormats.add("containsText", {
      text: "retained", format: { fill: "#DBEAFE", font: { color: "#1D4ED8", bold: true } },
    });
    sheet.getRange(`Q2:Q${rowCount}`).conditionalFormats.add("containsText", {
      text: "True", format: { fill: "#DCFCE7", font: { color: "#166534", bold: true } },
    });
    sheet.getRange(`Q2:Q${rowCount}`).conditionalFormats.add("containsText", {
      text: "False", format: { fill: "#FEE2E2", font: { color: "#991B1B" } },
    });
    const table = sheet.tables.add(
      sheet.getRangeByIndexes(0, 0, rowCount, lastColumn), true, tableName,
    );
    table.showFilterButton = true;
  }
}

styleDataSheet(all, lastRow, "AllEvaluationsTable");
styleDataSheet(current, currentRows.length + 1, "Current8MTable");

await fs.mkdir(previewDir, { recursive: true });
for (const [sheetName, range, fileName] of [
  ["Summary", "A1:H12", "summary.png"],
  ["8M 150t", `A1:N${Math.min(currentRows.length + 1, 25)}`, "current.png"],
  ["All Evaluations", "A1:N25", "all-evaluations.png"],
]) {
  const preview = await workbook.render({ sheetName, range, scale: 1.5, format: "png" });
  await fs.writeFile(`${previewDir}/${fileName}`, new Uint8Array(await preview.arrayBuffer()));
}

console.log((await workbook.inspect({
  kind: "table", range: "Summary!A1:H12", include: "values,formulas",
  tableMaxRows: 12, tableMaxCols: 8,
})).ndjson);
console.log((await workbook.inspect({
  kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 }, summary: "formula error scan",
})).ndjson);

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(`output=${outputPath}`);
