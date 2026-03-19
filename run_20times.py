#!/usr/bin/env python3
"""
Run exact (1 time) + V2 (20 times) on small test data, export comparison to Excel.
"""

import subprocess
import os
import re

# ============================================================
# CONFIG
# ============================================================

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

# Choose dataset: set DATASET env var or default to "tiny"
# Options: tiny, dataset_a, dataset_b, dataset_c, dataset_d
DATASET = os.environ.get("DATASET", "tiny")

RESULTS_DIR = os.path.join(PROJECT_DIR, "results", f"runs_{DATASET}")
NUM_RUNS = 20
NUM_TYPES = 18

NVERTS = os.path.join(PROJECT_DIR, "testdata", DATASET, "nverts.txt")
SIMPLICES = os.path.join(PROJECT_DIR, "testdata", DATASET, "simplices.txt")

# Fallback for old "tiny" format (flat testdata/)
if not os.path.exists(NVERTS):
    NVERTS = os.path.join(PROJECT_DIR, "testdata", "nverts.txt")
    SIMPLICES = os.path.join(PROJECT_DIR, "testdata", "simplices.txt")

# ============================================================
# SETUP
# ============================================================

os.makedirs(RESULTS_DIR, exist_ok=True)

# Build
subprocess.run(["make", "all"], cwd=PROJECT_DIR, capture_output=True)
subprocess.run(["make", "testdata"], cwd=PROJECT_DIR, capture_output=True)

# ============================================================
# PARSE .res FILE
# ============================================================

def parse_res(filepath):
    data = {"method": "", "time_ms": 0, "total": 0, "samples": 0,
            "counts": [0]*NUM_TYPES, "freq": [0.0]*NUM_TYPES,
            "ci_lo": [0.0]*NUM_TYPES, "ci_hi": [0.0]*NUM_TYPES}
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith("METHOD:"):
                data["method"] = line.split(":", 1)[1].strip()
            elif line.startswith("TIME_MS:"):
                data["time_ms"] = float(line.split(":", 1)[1].strip())
            elif line.startswith("TOTAL:"):
                data["total"] = int(line.split(":", 1)[1].strip())
            elif line.startswith("SAMPLES:"):
                data["samples"] = int(line.split(":", 1)[1].strip())
            elif line[0:1].isdigit():
                parts = line.split(",")
                if len(parts) >= 3:
                    t = int(parts[0]) - 1
                    if 0 <= t < NUM_TYPES:
                        data["counts"][t] = int(parts[1])
                        data["freq"][t] = float(parts[2])
                        if len(parts) >= 5:
                            data["ci_lo"][t] = float(parts[3])
                            data["ci_hi"][t] = float(parts[4])
    return data

# ============================================================
# RUN EXACT (once)
# ============================================================

print(f"Dataset: {DATASET}")
print(f"  nverts:    {NVERTS}")
print(f"  simplices: {SIMPLICES}")
print()

print("Running exact enumeration...")
exact_file = os.path.join(RESULTS_DIR, "exact.res")
subprocess.run(
    [os.path.join(PROJECT_DIR, "exact_sfd"), NVERTS, SIMPLICES, exact_file],
    cwd=PROJECT_DIR, capture_output=True
)
exact = parse_res(exact_file)
print(f"  Exact: {exact['total']} simplets, {exact['time_ms']:.3f} ms")

# ============================================================
# RUN V2 x 20
# ============================================================

v2_runs = []
for i in range(1, NUM_RUNS + 1):
    print(f"Running V2 #{i}/{NUM_RUNS}...", end=" ", flush=True)
    out_file = os.path.join(RESULTS_DIR, f"v2_run{i}.res")
    subprocess.run(
        [os.path.join(PROJECT_DIR, "unified_sfd_v2"), NVERTS, SIMPLICES, out_file, "0"],
        cwd=PROJECT_DIR, capture_output=True
    )
    data = parse_res(out_file)
    v2_runs.append(data)
    print(f"samples={data['samples']}, time={data['time_ms']:.0f}ms")

print(f"\nAll {NUM_RUNS} runs complete.")

# ============================================================
# BUILD EXCEL
# ============================================================

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side, numbers
from openpyxl.utils import get_column_letter

wb = Workbook()

# --- Colors ---
HEADER_FILL = PatternFill("solid", fgColor="2F5496")
HEADER_FONT = Font(name="Arial", bold=True, color="FFFFFF", size=11)
EXACT_FILL = PatternFill("solid", fgColor="D6E4F0")
V2_FILL = PatternFill("solid", fgColor="E2EFDA")
STAT_FILL = PatternFill("solid", fgColor="FFF2CC")
ERR_FILL = PatternFill("solid", fgColor="FCE4EC")
BOLD = Font(name="Arial", bold=True, size=11)
NORMAL = Font(name="Arial", size=11)
SMALL = Font(name="Arial", size=10)
CENTER = Alignment(horizontal="center", vertical="center")
THIN_BORDER = Border(
    left=Side(style="thin"), right=Side(style="thin"),
    top=Side(style="thin"), bottom=Side(style="thin")
)

def style_cell(cell, font=NORMAL, fill=None, align=CENTER, border=THIN_BORDER, fmt=None):
    cell.font = font
    if fill: cell.fill = fill
    cell.alignment = align
    cell.border = border
    if fmt: cell.number_format = fmt

def auto_width(ws):
    for col in ws.columns:
        max_len = 0
        col_letter = get_column_letter(col[0].column)
        for cell in col:
            if cell.value:
                max_len = max(max_len, len(str(cell.value)))
        ws.column_dimensions[col_letter].width = min(max(max_len + 3, 10), 20)

# ============================================================
# SHEET 1: Per-Run Frequencies
# ============================================================

ws1 = wb.active
ws1.title = "V2 Frequencies"

# Find active types (non-zero in exact or any v2 run)
active_types = []
for t in range(NUM_TYPES):
    if exact["freq"][t] > 0.0001 or any(r["freq"][t] > 0.0001 for r in v2_runs):
        active_types.append(t)

# Header row 1: Title
ws1.merge_cells(start_row=1, start_column=1, end_row=1, end_column=2 + len(active_types))
c = ws1.cell(row=1, column=1, value=f"SFD Comparison: Exact vs V2 (20 runs, N=0, {DATASET})")
style_cell(c, font=Font(name="Arial", bold=True, size=14, color="2F5496"))

# Header row 3
row = 3
style_cell(ws1.cell(row=row, column=1, value="Run"), HEADER_FONT, HEADER_FILL)
style_cell(ws1.cell(row=row, column=2, value="Samples"), HEADER_FONT, HEADER_FILL)
for ci, t in enumerate(active_types):
    style_cell(ws1.cell(row=row, column=3+ci, value=f"Type {t+1}"), HEADER_FONT, HEADER_FILL)
style_cell(ws1.cell(row=row, column=3+len(active_types), value="Time (ms)"), HEADER_FONT, HEADER_FILL)

# Exact row
row = 4
style_cell(ws1.cell(row=row, column=1, value="EXACT"), BOLD, EXACT_FILL)
style_cell(ws1.cell(row=row, column=2, value=exact["total"]), NORMAL, EXACT_FILL)
for ci, t in enumerate(active_types):
    style_cell(ws1.cell(row=row, column=3+ci, value=exact["freq"][t]), NORMAL, EXACT_FILL, fmt="0.000000")
style_cell(ws1.cell(row=row, column=3+len(active_types), value=exact["time_ms"]), NORMAL, EXACT_FILL, fmt="0.000")

# V2 runs
for i, run in enumerate(v2_runs):
    row = 5 + i
    style_cell(ws1.cell(row=row, column=1, value=f"V2 #{i+1}"), NORMAL, V2_FILL)
    style_cell(ws1.cell(row=row, column=2, value=run["samples"]), NORMAL, V2_FILL)
    for ci, t in enumerate(active_types):
        style_cell(ws1.cell(row=row, column=3+ci, value=run["freq"][t]), NORMAL, fmt="0.000000")
    style_cell(ws1.cell(row=row, column=3+len(active_types), value=run["time_ms"]), NORMAL, fmt="0.000")

# Statistics rows
stat_row = 5 + NUM_RUNS + 1
labels = ["MEAN", "STD", "MIN", "MAX", "EXACT", "MEAN Error"]
for si, label in enumerate(labels):
    r = stat_row + si
    style_cell(ws1.cell(row=r, column=1, value=label), BOLD, STAT_FILL)
    style_cell(ws1.cell(row=r, column=2, value=""), NORMAL, STAT_FILL)
    
    first_data = 5  # row of V2 #1
    last_data = 5 + NUM_RUNS - 1  # row of V2 #20
    
    for ci, t in enumerate(active_types):
        col_letter = get_column_letter(3 + ci)
        c = ws1.cell(row=r, column=3+ci)
        style_cell(c, NORMAL, STAT_FILL, fmt="0.000000")
        
        if label == "MEAN":
            c.value = f"=AVERAGE({col_letter}{first_data}:{col_letter}{last_data})"
        elif label == "STD":
            c.value = f"=STDEV({col_letter}{first_data}:{col_letter}{last_data})"
        elif label == "MIN":
            c.value = f"=MIN({col_letter}{first_data}:{col_letter}{last_data})"
        elif label == "MAX":
            c.value = f"=MAX({col_letter}{first_data}:{col_letter}{last_data})"
        elif label == "EXACT":
            c.value = exact["freq"][t]
        elif label == "MEAN Error":
            exact_row = stat_row + 4  # EXACT stats row
            mean_row = stat_row       # MEAN stats row
            c.value = f"=ABS({col_letter}{mean_row}-{col_letter}{exact_row})"
            style_cell(c, NORMAL, ERR_FILL, fmt="0.000000")
    
    # Time column stats
    time_col = get_column_letter(3 + len(active_types))
    tc = ws1.cell(row=r, column=3+len(active_types))
    style_cell(tc, NORMAL, STAT_FILL, fmt="0.000")
    if label == "MEAN":
        tc.value = f"=AVERAGE({time_col}{first_data}:{time_col}{last_data})"
    elif label == "STD":
        tc.value = f"=STDEV({time_col}{first_data}:{time_col}{last_data})"
    elif label == "MIN":
        tc.value = f"=MIN({time_col}{first_data}:{time_col}{last_data})"
    elif label == "MAX":
        tc.value = f"=MAX({time_col}{first_data}:{time_col}{last_data})"
    elif label == "EXACT":
        tc.value = exact["time_ms"]
    elif label == "MEAN Error":
        tc.value = ""

auto_width(ws1)

# ============================================================
# SHEET 2: Per-Run Absolute Errors
# ============================================================

ws2 = wb.create_sheet("V2 Errors")

ws2.merge_cells(start_row=1, start_column=1, end_row=1, end_column=2 + len(active_types))
c = ws2.cell(row=1, column=1, value="Absolute Error: |V2 - Exact| per run")
style_cell(c, font=Font(name="Arial", bold=True, size=14, color="2F5496"))

row = 3
style_cell(ws2.cell(row=row, column=1, value="Run"), HEADER_FONT, HEADER_FILL)
for ci, t in enumerate(active_types):
    style_cell(ws2.cell(row=row, column=2+ci, value=f"Type {t+1}"), HEADER_FONT, HEADER_FILL)
style_cell(ws2.cell(row=row, column=2+len(active_types), value="RMSE"), HEADER_FONT, HEADER_FILL)
style_cell(ws2.cell(row=row, column=3+len(active_types), value="MAE"), HEADER_FONT, HEADER_FILL)

for i, run in enumerate(v2_runs):
    row = 4 + i
    style_cell(ws2.cell(row=row, column=1, value=f"V2 #{i+1}"), NORMAL, V2_FILL)
    
    errors = []
    for ci, t in enumerate(active_types):
        err = abs(run["freq"][t] - exact["freq"][t])
        errors.append(err)
        c = ws2.cell(row=row, column=2+ci, value=err)
        fill = ERR_FILL if err > 0.01 else None
        style_cell(c, NORMAL, fill, fmt="0.000000")
    
    # RMSE
    rmse = (sum(e**2 for e in errors) / len(errors)) ** 0.5
    style_cell(ws2.cell(row=row, column=2+len(active_types), value=rmse), NORMAL, fmt="0.000000")
    
    # MAE
    mae = sum(errors) / len(errors)
    style_cell(ws2.cell(row=row, column=3+len(active_types), value=mae), NORMAL, fmt="0.000000")

# Stats
stat_row = 4 + NUM_RUNS + 1
for si, label in enumerate(["MEAN", "STD", "MAX"]):
    r = stat_row + si
    style_cell(ws2.cell(row=r, column=1, value=label), BOLD, STAT_FILL)
    
    first = 4
    last = 4 + NUM_RUNS - 1
    
    for ci in range(len(active_types) + 2):  # +2 for RMSE, MAE
        col_letter = get_column_letter(2 + ci)
        c = ws2.cell(row=r, column=2+ci)
        style_cell(c, NORMAL, STAT_FILL, fmt="0.000000")
        if label == "MEAN":
            c.value = f"=AVERAGE({col_letter}{first}:{col_letter}{last})"
        elif label == "STD":
            c.value = f"=STDEV({col_letter}{first}:{col_letter}{last})"
        elif label == "MAX":
            c.value = f"=MAX({col_letter}{first}:{col_letter}{last})"

auto_width(ws2)

# ============================================================
# SHEET 3: Confidence Interval Coverage
# ============================================================

ws3 = wb.create_sheet("CI Coverage")

ws3.merge_cells(start_row=1, start_column=1, end_row=1, end_column=6)
c = ws3.cell(row=1, column=1, value="95% CI Coverage: Does CI contain exact value?")
style_cell(c, font=Font(name="Arial", bold=True, size=14, color="2F5496"))

row = 3
for h, val in enumerate(["Run", "Type", "Exact", "CI Low", "CI High", "Covered?"]):
    style_cell(ws3.cell(row=row, column=1+h, value=val), HEADER_FONT, HEADER_FILL)

row = 4
total_checks = 0
total_covered = 0

for i, run in enumerate(v2_runs):
    for t in active_types:
        ex = exact["freq"][t]
        lo = run["ci_lo"][t]
        hi = run["ci_hi"][t]
        covered = lo <= ex <= hi
        total_checks += 1
        if covered: total_covered += 1
        
        style_cell(ws3.cell(row=row, column=1, value=f"V2 #{i+1}"), SMALL)
        style_cell(ws3.cell(row=row, column=2, value=f"Type {t+1}"), SMALL)
        style_cell(ws3.cell(row=row, column=3, value=ex), SMALL, fmt="0.000000")
        style_cell(ws3.cell(row=row, column=4, value=lo), SMALL, fmt="0.000000")
        style_cell(ws3.cell(row=row, column=5, value=hi), SMALL, fmt="0.000000")
        
        c = ws3.cell(row=row, column=6, value="YES" if covered else "NO")
        fill = PatternFill("solid", fgColor="C6EFCE") if covered else PatternFill("solid", fgColor="FFC7CE")
        style_cell(c, Font(name="Arial", bold=True, size=10, color="006100" if covered else "9C0006"), fill)
        row += 1

# Summary
row += 1
ws3.merge_cells(start_row=row, start_column=1, end_row=row, end_column=6)
pct = 100.0 * total_covered / total_checks if total_checks > 0 else 0
c = ws3.cell(row=row, column=1,
    value=f"Overall Coverage: {total_covered}/{total_checks} = {pct:.1f}%")
style_cell(c, Font(name="Arial", bold=True, size=13, color="2F5496"))

auto_width(ws3)

# ============================================================
# SHEET 4: Summary
# ============================================================

ws4 = wb.create_sheet("Summary")
ws4.sheet_properties.tabColor = "2F5496"

# Move Summary to first position
wb.move_sheet("Summary", offset=-3)

ws4.merge_cells("A1:D1")
style_cell(ws4.cell(row=1, column=1, value=f"SFD Estimation: Exact vs V2 ({DATASET})"),
           Font(name="Arial", bold=True, size=16, color="2F5496"))

row = 3
summary_data = [
    ("Dataset", DATASET),
    ("N (subgraph)", "0 (use all)"),
    ("V2 Runs", str(NUM_RUNS)),
    ("V2 Mode", "Paper mode (size 2-4 only)"),
    ("Exact Total Simplets", str(exact["total"])),
    ("", ""),
]

for label, val in summary_data:
    style_cell(ws4.cell(row=row, column=1, value=label), BOLD)
    style_cell(ws4.cell(row=row, column=2, value=val), NORMAL)
    row += 1

# Per-type summary table
row += 1
style_cell(ws4.cell(row=row, column=1, value="Type"), HEADER_FONT, HEADER_FILL)
style_cell(ws4.cell(row=row, column=2, value="Exact"), HEADER_FONT, HEADER_FILL)
style_cell(ws4.cell(row=row, column=3, value="V2 Mean"), HEADER_FONT, HEADER_FILL)
style_cell(ws4.cell(row=row, column=4, value="V2 Std"), HEADER_FONT, HEADER_FILL)
style_cell(ws4.cell(row=row, column=5, value="Mean Error"), HEADER_FONT, HEADER_FILL)
style_cell(ws4.cell(row=row, column=6, value="CI Coverage"), HEADER_FONT, HEADER_FILL)

for t in active_types:
    row += 1
    freqs = [r["freq"][t] for r in v2_runs]
    mean_f = sum(freqs) / len(freqs)
    std_f = (sum((f - mean_f)**2 for f in freqs) / len(freqs)) ** 0.5
    mean_err = abs(mean_f - exact["freq"][t])
    
    ci_covered = sum(1 for r in v2_runs if r["ci_lo"][t] <= exact["freq"][t] <= r["ci_hi"][t])
    ci_pct = 100.0 * ci_covered / NUM_RUNS
    
    style_cell(ws4.cell(row=row, column=1, value=f"Type {t+1}"), BOLD)
    style_cell(ws4.cell(row=row, column=2, value=exact["freq"][t]), NORMAL, EXACT_FILL, fmt="0.000000")
    style_cell(ws4.cell(row=row, column=3, value=mean_f), NORMAL, V2_FILL, fmt="0.000000")
    style_cell(ws4.cell(row=row, column=4, value=std_f), NORMAL, fmt="0.000000")
    style_cell(ws4.cell(row=row, column=5, value=mean_err), NORMAL, ERR_FILL if mean_err > 0.01 else None, fmt="0.000000")
    style_cell(ws4.cell(row=row, column=6, value=f"{ci_pct:.0f}%"), NORMAL,
               PatternFill("solid", fgColor="C6EFCE") if ci_pct >= 90 else ERR_FILL)

# Overall stats
row += 2
all_mean_errors = []
for t in active_types:
    freqs = [r["freq"][t] for r in v2_runs]
    mean_f = sum(freqs) / len(freqs)
    all_mean_errors.append(abs(mean_f - exact["freq"][t]))

overall_rmse = (sum(e**2 for e in all_mean_errors) / len(all_mean_errors)) ** 0.5
overall_mae = sum(all_mean_errors) / len(all_mean_errors)
overall_ci = 100.0 * total_covered / total_checks if total_checks > 0 else 0

style_cell(ws4.cell(row=row, column=1, value="Overall RMSE (mean)"), BOLD)
style_cell(ws4.cell(row=row, column=2, value=overall_rmse), NORMAL, fmt="0.000000")
row += 1
style_cell(ws4.cell(row=row, column=1, value="Overall MAE (mean)"), BOLD)
style_cell(ws4.cell(row=row, column=2, value=overall_mae), NORMAL, fmt="0.000000")
row += 1
style_cell(ws4.cell(row=row, column=1, value="Overall CI Coverage"), BOLD)
style_cell(ws4.cell(row=row, column=2, value=f"{overall_ci:.1f}%"), NORMAL,
           PatternFill("solid", fgColor="C6EFCE") if overall_ci >= 90 else ERR_FILL)
row += 1
avg_time = sum(r["time_ms"] for r in v2_runs) / NUM_RUNS
style_cell(ws4.cell(row=row, column=1, value="Avg V2 Time (ms)"), BOLD)
style_cell(ws4.cell(row=row, column=2, value=avg_time), NORMAL, fmt="0.0")

auto_width(ws4)

# ============================================================
# SAVE
# ============================================================

output_path = os.path.join(PROJECT_DIR, "results", f"sfd_exact_vs_v2_{DATASET}_20runs.xlsx")
os.makedirs(os.path.dirname(output_path), exist_ok=True)
wb.save(output_path)
print(f"\nExcel saved: {output_path}")