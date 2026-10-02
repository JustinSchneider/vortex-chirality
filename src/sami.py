"""SAMI DR3 access for the confirmatory test C1 (preregistration/OSF_C1_confirmatory.md).

Catalogues (non-kinematic columns only) were fetched from the Data Central
TAP service into data/sami/*.csv. Maps are the Data Central bulk-download
products, unpacked under data/sami/maps/ (any directory layout):
  stars: two-moment pPXF, default 0.5" spaxels, velocity (VEL, VEL_ERR, SN)
         and dispersion (SIG, SIG_ERR, SN);
  gas:   LZIFU 1-comp velocity (+ V_ERR), dispersion (+ VDISP_ERR) and
         H-alpha flux (+ error) maps.
DR3 file names are not documented, so files are located by keyword
(FILE_KEYS) and the match must be unique; the dry run checks this.

Axes follow the MaNGA pipeline: x along increasing array column, y along
increasing row, origin at the map centre (SAMI maps are centred on the
galaxy; galaxies flagged WARNWCS are excluded). Nothing here prints or
summarises rotation amplitudes.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from src.manga import inclination_deg

ROOT = Path(__file__).resolve().parents[1]
SAMI = ROOT / "data" / "sami"
MAPS = SAMI / "maps"
SPAXEL_ARCSEC = 0.5
MAX_VERR = 30.0
MIN_HA_SNR = 5.0
GOOD_BAD_CLASS = (0, 5, 8)            # Croom et al. 2021, Sect. 2.4
COSMO_H0 = 70.0

# keyword sets identifying each product in a file name (all must appear;
# none of the 'not' words may appear).
FILE_KEYS = {
    "stellar_vel": (("stellar", "velocity", "default", "two-moment"), ("dispersion",)),
    "stellar_sig": (("stellar", "dispersion", "default", "two-moment"), ()),
    "gas_vel": (("gas", "velocity", "default", "1-comp"), ("dispersion", "vdisp")),
    "gas_sig": (("gas", "vdisp", "default", "1-comp"), ()),
    "ha_flux": (("halpha", "default", "1-comp"), ()),
}


def catalogue() -> pd.DataFrame:
    """One row per SAMI DR3 galaxy with a best cube."""
    g = pd.read_csv(SAMI / "gama.csv")
    c = pd.read_csv(SAMI / "clusters.csv")
    g["z"] = g["z_tonry"].where(g["z_tonry"] > 0, g["z_spec"])
    c["z"] = c["z_spec"]
    cols = ["CATID", "RA_OBJ", "DEC_OBJ", "z", "r_e", "ellip", "PA", "Mstar", "BAD_CLASS"]
    d = pd.concat([g[cols], c[cols]]).drop_duplicates("CATID")
    cube = pd.read_csv(SAMI / "cubeobs.csv")
    cube = cube[cube["ISBEST"] == 1].drop_duplicates("CATID")
    d = d.merge(cube, on="CATID", how="inner")
    m = pd.read_csv(SAMI / "morph.csv")
    d = d.merge(m, on="CATID", how="left")
    d = d.rename(columns={"RA_OBJ": "ra", "DEC_OBJ": "dec", "r_e": "Re_arcsec",
                          "Mstar": "logMstar", "PA": "pa_phot"})
    d["morph"] = d["TYPE"].where(d["TYPE"].between(0, 3))
    d["ba"] = 1 - d["ellip"]
    d["inc"] = inclination_deg(d["ba"])
    # angular-diameter distance, flat LCDM (Om = 0.3, H0 = 70)
    from astropy.cosmology import FlatLambdaCDM
    da = FlatLambdaCDM(H0=COSMO_H0, Om0=0.3).angular_diameter_distance(
        d["z"].clip(lower=1e-4).values).value
    d["kpc_per_arcsec"] = da * 1e3 * np.pi / (180 * 3600)
    d["survey"] = "SAMI"
    d["quality_ok"] = (d["BAD_CLASS"].isin(GOOD_BAD_CLASS) & (d["WARNSTAR"] == 0)
                       & (d["WARNSK2M"] == 0) & (d["WARNSKER"] == 0)
                       & (d["WARNWCS"] == 0) & (d["WARNEMFT"] == 0)
                       & (d["WARNMULT"] == 0))
    d["name"] = d["CATID"].astype(str)
    return d.set_index("name")


def find_file(catid, product, files=None):
    """Unique file for (CATID, product) under data/sami/maps."""
    want, avoid = FILE_KEYS[product]
    files = files if files is not None else list(MAPS.rglob(f"{catid}_*.fits*"))
    hits = [f for f in files
            if f.name.split("_")[0] == str(catid)
            and all(w in f.name.lower() for w in want)
            and not any(a in f.name.lower() for a in avoid)]
    if len(hits) != 1:
        raise FileNotFoundError(f"{catid} {product}: {len(hits)} matches")
    return hits[0]


def _read(path, ext=0):
    with fits.open(path) as h:
        return np.asarray(h[ext].data, float)


def _ext(path, name):
    with fits.open(path) as h:
        return np.asarray(h[name].data, float)


def available(catid):
    try:
        for p in FILE_KEYS:
            find_file(catid, p)
        return True
    except FileNotFoundError:
        return False


def galaxy_arrays(catid, z=None, session=None):
    """Masked stellar and gas maps; format of g1_manga_test.galaxy_arrays."""
    f = {p: find_file(catid, p) for p in FILE_KEYS}
    sv, sve, ssn = _read(f["stellar_vel"]), _ext(f["stellar_vel"], "VEL_ERR"), _ext(f["stellar_vel"], "SN")
    ss, sse = _read(f["stellar_sig"]), _ext(f["stellar_sig"], "SIG_ERR")
    gv, gve = _read(f["gas_vel"]), _ext(f["gas_vel"], "V_ERR")
    gs = _read(f["gas_sig"])
    ha, hae = _read(f["ha_flux"]), _ext(f["ha_flux"], 1)
    # 1-comp maps are 2-D; tolerate a leading length-1 axis
    gv, gve, gs, ha, hae = (a[0] if a.ndim == 3 else a for a in (gv, gve, gs, ha, hae))
    ny, nx = sv.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    x = (xx - (nx - 1) / 2) * SPAXEL_ARCSEC
    y = (yy - (ny - 1) / 2) * SPAXEL_ARCSEC
    with np.errstate(invalid="ignore", divide="ignore"):
        s_ok = np.isfinite(sv) & (sve > 0) & (sve < MAX_VERR) & (ssn > 3)
        sig_ok = s_ok & np.isfinite(ss) & (sse < 0.1 * ss + 25) & (ss > 35)
        snr = np.where(hae > 0, ha / hae, 0.0)
        g_ok = np.isfinite(gv) & (gve > 0) & (gve <= MAX_VERR) & (snr >= MIN_HA_SNR)
    nan = np.nan
    gsig = np.where(np.isfinite(gs) & (gs > 0), gs, 0.0)   # LZIFU: instrument-corrected
    maps = dict(
        vs=np.where(s_ok, sv - np.nanmedian(sv[s_ok]) if s_ok.any() else nan, nan),
        evs=np.where(s_ok, sve, np.inf),
        ss=np.where(sig_ok, ss, nan),
        vg=np.where(g_ok, gv - np.nanmedian(gv[g_ok]) if g_ok.any() else nan, nan),
        evg=np.where(g_ok, gve, np.inf),
        sg=np.where(g_ok, gsig, nan))
    return maps, x, y, s_ok, g_ok
