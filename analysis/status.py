#!/usr/bin/env python
"""Status of the MaNGA G1 pipeline: what is running, and how far along.

Usage (from C:\\Science\\VORTEX):  python analysis/status.py

Stages: 1 selection (g1_manga_select.py) -> 2 registered test
(g1_manga_test.py) -> 3 exploratory inflow (g1_manga_inflow_exploratory.py).
Selection progress shows counts only, never rotation speeds.
"""
import re
import subprocess
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "data" / "processed" / "manga_pa.csv"
TOTAL = 7582
STAGES = [
    ("1 selection", "g1_manga_select.py", ROOT / "results" / "g1_manga_select.log"),
    ("2 registered test", "g1_manga_test.py", ROOT / "results" / "g1_manga_test.log"),
    ("3 inflow (exploratory)", "g1_manga_inflow_exploratory.py",
     ROOT / "results" / "g1_manga_inflow.log"),
]


def running_scripts():
    """Names of pipeline scripts with a live python process (Windows)."""
    cmd = ("Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
           "ForEach-Object { $_.CommandLine }")
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-Command", cmd],
                             capture_output=True, text=True, timeout=60).stdout
    except Exception:
        return None
    return out


def ago(path):
    if not path.exists():
        return "never"
    m = (time.time() - path.stat().st_mtime) / 60
    return f"{m:.0f} min ago" if m < 120 else f"{m / 60:.1f} h ago"


procs = running_scripts()
print("=" * 60)
print("MaNGA G1 pipeline status   " + time.strftime("%Y-%m-%d %H:%M"))
print("=" * 60)
for name, script, log in STAGES:
    alive = procs is not None and script in procs
    n_workers = procs.count("multiprocessing-fork") if alive and procs else 0
    state = "RUNNING" if alive else "not running"
    extra = f" ({n_workers} worker processes)" if alive and n_workers else ""
    print(f"Stage {name:24s}: {state}{extra}; log last written {ago(log)}")

# Stage 1 detail
if CSV.exists():
    d = pd.read_csv(CSV).drop_duplicates("plateifu", keep="last")
    ok = d[(d["e_pa_star"] <= 20) & (d["e_pa_gas"] <= 20)]
    print(f"\n[1] galaxies processed: {len(d)} / {TOTAL} ({100 * len(d) / TOTAL:.1f}%)")
    print(f"    good PAs: {len(ok)}; counter-rotator candidates:"
          f" {(ok['dpa'] > 150).sum()}; co-rotator candidates: {(ok['dpa'] < 30).sum()}")
    sel_log = STAGES[0][2]
    if sel_log.exists():
        m = re.findall(r"(\d+)/(\d+) processed", sel_log.read_text(errors="replace"))
        if m:
            print(f"    last progress line: {m[-1][0]}/{m[-1][1]}")

# Stage 2/3 detail: show the log as it stands (written as each step finishes)
for name, script, log in STAGES[1:]:
    if log.exists() and (not CSV.exists() or log.stat().st_mtime >= CSV.stat().st_mtime):
        txt = log.read_text(errors="replace").strip()
        print(f"\n[{name}] log so far:")
        print("    " + ("\n    ".join(txt.splitlines()) if txt else "(started; nothing written yet)"))
    elif log.exists():
        print(f"\n[{name}] log on disk is from an EARLIER run (older than the selection);"
              " the new run has not written yet.")

if procs is not None and not any(s in procs for _, s, _ in STAGES):
    print("\nNo pipeline process is running. If the last stage's log ends with a"
          " verdict/result, the chain finished; otherwise it stopped early.")
