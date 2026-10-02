#!/usr/bin/env python
"""Download CALIFA products for the C1 parent sample (no measurement, no output
beyond progress counts). Resumable."""
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis.c1_test import parent_sample  # noqa: E402
from src import califa  # noqa: E402


def one(name):
    califa.ensure_cached(name, requests.Session())
    return name


if __name__ == "__main__":
    d, _ = parent_sample()
    names = [n for n in d[d.survey == "CALIFA"].index
             if not (califa.CACHE / f"{n}_gas.npz").exists()
             or not (califa.CACHE / f"{n}.CALIFA.V1200.stekin.fits").exists()]
    print(f"CALIFA to fetch: {len(names)}", flush=True)
    fails = 0
    with ThreadPoolExecutor(3) as ex:
        futs = {ex.submit(one, n): n for n in names}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                f.result()
            except Exception as e:  # noqa: BLE001
                fails += 1
                print(f"  fail {futs[f]}: {type(e).__name__}", flush=True)
            if i % 10 == 0 or i == len(names):
                print(f"  {i}/{len(names)} ({fails} failed)", flush=True)
