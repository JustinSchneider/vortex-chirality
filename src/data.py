"""Loaders for SPARC and the existing Rational Taper fit tables.

All files live in data/raw (copied unchanged from
RT/rational-taper-validation/data). Nothing here modifies the source data.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

# Mass-to-light ratios at 3.6 micron (Lelli+2016; same as RT pipeline)
UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7
HELIUM = 1.33

# SPARC Table 1 columns. The published byte ranges are offset by one from
# the actual rows, but every field is present and space-free, so the rows
# are parsed by whitespace splitting.
_SPARC_COLS = [
    ("Galaxy", str), ("T", int), ("D", float), ("e_D", float),
    ("f_D", int), ("Inc", float), ("e_Inc", float), ("L36", float),
    ("e_L36", float), ("Reff", float), ("SBeff", float), ("Rdisk", float),
    ("SBdisk", float), ("MHI", float), ("RHI", float), ("Vflat", float),
    ("e_Vflat", float), ("Q", int), ("Ref", str),
]


def load_sparc_table() -> pd.DataFrame:
    """SPARC Table 1 (Lelli+2016c), one row per galaxy."""
    lines = (RAW / "SPARC_Lelli2016c.mrt").read_text().splitlines()
    start = max(i for i, l in enumerate(lines) if l.startswith("----")) + 1
    rows = []
    for line in lines[start:]:
        tok = line.split()
        if len(tok) != len(_SPARC_COLS):
            continue
        rows.append({name: typ(v) for (name, typ), v in zip(_SPARC_COLS, tok)})
    df = pd.DataFrame(rows).set_index("Galaxy")
    if len(df) != 175:
        raise ValueError(f"Expected 175 SPARC galaxies, parsed {len(df)}")
    return df


def load_bulges() -> pd.Series:
    """Bulge luminosity (1e9 Lsun) per galaxy, from SPARC Bulges.mrt."""
    out = {}
    for line in (RAW / "Bulges.mrt").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, val = line.split()
        out[name] = float(val)
    return pd.Series(out, name="Lbul")


def load_mass_models() -> pd.DataFrame:
    """SPARC per-radius mass models (Lelli+2016c)."""
    from astropy.io import ascii
    t = ascii.read(RAW / "MassModels_Lelli2016c.mrt", format="mrt")
    return t.to_pandas()


def v_bary(vgas, vdisk, vbul, ud=UPSILON_DISK, ub=UPSILON_BULGE):
    """Signed-square baryonic speed (Lelli+2016 Eq. 2)."""
    v2 = (np.abs(vgas) * vgas + ud * np.abs(vdisk) * vdisk
          + ub * np.abs(vbul) * vbul)
    return np.sign(v2) * np.sqrt(np.abs(v2))


def load_rt_fits() -> pd.DataFrame:
    """Schneider 2026a fit table (Linear and Tapered fits, 175 galaxies)."""
    return pd.read_csv(RAW / "Schneider_2026_SPARC_Fit_Parameters.csv",
                       index_col="GalaxyID")


def load_tournament() -> pd.DataFrame:
    """Schneider 2026b model tournament (RT, NFW, MOND fits)."""
    return pd.read_csv(RAW / "Tournament_Results.csv", index_col="GalaxyID")
