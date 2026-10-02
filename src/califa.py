"""CALIFA access for the confirmatory test C1 (preregistration/OSF_C1_confirmatory.md).

Stars: Falcon-Barroso et al. (2017) V1200 Voronoi-binned pPXF tables, one
row per 1" spaxel (X, Y in arcsec relative to the galaxy centre).
Gas: eCALIFA pyPipe3D V500 v2.3 cubes (Sanchez et al. 2023), 0.5" spaxels
with a celestial WCS. H-alpha velocity/flux and errors from FLUX_ELINES;
the H-alpha line width (FWHM in Angstrom, instrumental width included) from
ELINES channel 1.

Both are placed on the gas 0.5" grid. Axes follow the MaNGA pipeline: x
along increasing array column (west), y along increasing row (north),
origin at the catalogue galaxy position. Each grid spaxel takes the stellar
values of the nearest 1" stellar spaxel within 0.75".

Catalogue access reads only named, non-kinematic columns. Nothing here
prints or summarises rotation amplitudes.
"""

import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree

from src.manga import inclination_deg

ROOT = Path(__file__).resolve().parents[1]
CALIFA = ROOT / "data" / "califa"
CACHE = CALIFA / "maps"
STEKIN_URL = ("https://califa.caha.es/FTP-PUB/dataproducts/Stellar_Kinematics_V1200/"
              "{name}.CALIFA.V1200.stekin.fits")
PIPE3D_URL = "https://ifs.astroscu.unam.mx/CALIFA/V500/v2.3/pyPipe3D/{name}.Pipe3D.cube.fits.gz"
SPAXEL_ARCSEC = 0.5
C_KMS = 299792.458
HA_REST = 6562.8
SIG_INST_AA = 2.6          # V500 instrumental sigma (6.0 A FWHM / 2.354)
MAX_VERR = 30.0
MIN_HA_SNR = 5.0
STAR_MATCH_ARCSEC = 0.75
# FLUX_ELINES plane names (header NAMEnnn) for H-alpha.
HA_NAMES = {"flux": "flux Ha", "vel": "vel Ha", "eflux": "e_flux Ha", "evel": "e_vel Ha"}

# Hubble type -> common 0-3 morphology scale (PREREG C1, "Measured variables").
MORPH = {**{f"E{i}": 0.0 for i in range(8)}, "S0": 1.0, "S0a": 1.0,
         "Sa": 2.0, "Sab": 2.0, "Sb": 2.0,
         "Sbc": 3.0, "Sc": 3.0, "Scd": 3.0, "Sd": 3.0, "Sdm": 3.0, "Sm": 3.0, "I": 3.0}


def stekin_names():
    return sorted(l.split(".CALIFA")[0] for l in
                  (CALIFA / "stekin_list.txt").read_text().split())


def catalogue() -> pd.DataFrame:
    """One row per eCALIFA galaxy with a V1200 stellar-kinematics table."""
    g = fits.getdata(CALIFA / "galaxies_properties.fits.gz", 1)
    a = pd.DataFrame({
        "name": pd.Series(g["cubename"].astype(str)).str.strip(),
        "ra": g["RA"].astype(float), "dec": g["DEC"].astype(float),
        "z": g["z"].astype(float),
        "hubble": pd.Series(g["type"].astype(str)).str.strip(),
        "qc_flag": np.asarray(g["QC_flag"]).astype(int),
    })
    p = fits.getdata(CALIFA / "eCALIFA.pyPipe3D.fits", 1)
    b = pd.DataFrame({
        "name": pd.Series(p["cubename"].astype(str)).str.strip(),
        "Re_arcsec": p["Re_arc"].astype(float), "ellip": p["ellip"].astype(float),
        "pa_phot": p["PA"].astype(float),
        # 'DA' is the angular scale in kpc/arcsec (= Re_kpc / Re_arc for every row)
        "kpc_per_arcsec": p["DA"].astype(float),
        "logMstar": p["log_Mass"].astype(float),
    })
    d = a.merge(b, on="name", how="inner")
    d = d[d["name"].isin(set(stekin_names()))].drop_duplicates("name")
    d["morph"] = d["hubble"].map(MORPH)
    # 'ellip' is an eccentricity (tracks 'ecc' in galaxies_properties, r = 0.998;
    # ellipticals rounder than spirals only under this reading): b/a = sqrt(1 - e^2)
    d["ba"] = np.sqrt(1 - d["ellip"].clip(0, 1) ** 2)
    d["inc"] = inclination_deg(d["ba"])
    d["survey"] = "CALIFA"
    d["quality_ok"] = d["qc_flag"] != 1
    return d.set_index("name")


def _download(url, dest, session, timeout=600):
    for attempt in range(5):
        try:
            with session.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            return
        except Exception:
            if attempt == 4:
                raise
            time.sleep(10 * (attempt + 1))


def _plane(hdr, label):
    for k in hdr:
        if k.startswith("NAME") and str(hdr[k]).strip() == label:
            return int(k[4:])
    raise KeyError(label)


def _extract_pipe3d(path, out):
    with fits.open(path) as h:
        fe, fh = h["FLUX_ELINES"].data, h["FLUX_ELINES"].header
        planes = {k: np.asarray(fe[_plane(fh, v)], float) for k, v in HA_NAMES.items()}
        wcs = WCS(h["ORG_HDR"].header).celestial
        np.savez_compressed(
            out, fwhm=np.asarray(h["ELINES"].data[1], float),
            gaia_mask=np.asarray(h["GAIA_MASK"].data, float),
            wcs_header=np.array(wcs.to_header_string()), **planes)


def ensure_cached(name, session=None):
    """Fetch the stellar table and the H-alpha planes of the pyPipe3D cube."""
    CACHE.mkdir(parents=True, exist_ok=True)
    s = session or requests.Session()
    st = CACHE / f"{name}.CALIFA.V1200.stekin.fits"
    if not st.exists():
        _download(STEKIN_URL.format(name=name), st, s)
    gz = CACHE / f"{name}_gas.npz"
    if not gz.exists():
        local = CACHE / f"{name}.Pipe3D.cube.fits.gz"
        if local.exists():
            _extract_pipe3d(local, gz)
        else:
            with tempfile.TemporaryDirectory() as td:
                f = Path(td) / "cube.fits.gz"
                _download(PIPE3D_URL.format(name=name), f, s)
                _extract_pipe3d(f, gz)
    return st, gz


def galaxy_arrays(name, ra, dec, z, session=None):
    """Masked stellar and gas maps on the common 0.5" grid.

    Returns (maps, x, y, s_ok, g_ok) in the format of
    analysis.g1_manga_test.galaxy_arrays.
    """
    st_path, gz_path = ensure_cached(name, session)
    g = np.load(gz_path)
    wcs = WCS(fits.Header.fromstring(str(g["wcs_header"])))
    col0, row0 = wcs.world_to_pixel_values(ra, dec)
    ny, nx = g["vel"].shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    x = (xx - col0) * SPAXEL_ARCSEC
    y = (yy - row0) * SPAXEL_ARCSEC

    # gas
    ev, f, ef = g["evel"], g["flux"], g["eflux"]
    with np.errstate(invalid="ignore", divide="ignore"):
        snr = np.where(ef > 0, f / ef, 0.0)
    g_ok = (np.isfinite(g["vel"]) & (g["vel"] != 0) & (ev > 0) & (ev <= MAX_VERR)
            & (snr >= MIN_HA_SNR) & (g["gaia_mask"] == 0))
    lam = HA_REST * (1 + z)
    s2 = (g["fwhm"] / 2.354 / lam) ** 2 - (SIG_INST_AA / lam) ** 2
    g_sig = np.where(s2 > 0, C_KMS * np.sqrt(np.abs(s2)), 0.0)   # below resolution: cold

    # stars: nearest 1" spaxel
    t = fits.getdata(st_path, 1)
    X, Y = np.asarray(t["X"], float), np.asarray(t["Y"], float)
    tree = cKDTree(np.column_stack([X, Y]))
    dist, idx = tree.query(np.column_stack([x.ravel(), y.ravel()]))
    near = (dist <= STAR_MATCH_ARCSEC).reshape(x.shape)
    idx = idx.reshape(x.shape)

    def star(col):
        return np.where(near, np.asarray(t[col], float)[idx], np.nan)

    vp, dvp, sp, dsp = star("Vp"), star("DVp"), star("Sp"), star("DSp")
    s_ok = near & np.isfinite(vp) & (dvp > 0) & (dvp <= MAX_VERR)
    sig_ok = s_ok & (dsp < 0.1 * sp + 25) & (sp > 35)

    nan = np.nan
    maps = dict(
        vs=np.where(s_ok, vp - np.nanmedian(vp[s_ok]) if s_ok.any() else nan, nan),
        evs=np.where(s_ok, dvp, np.inf),
        ss=np.where(sig_ok, sp, nan),
        vg=np.where(g_ok, g["vel"] - np.nanmedian(g["vel"][g_ok]) if g_ok.any() else nan, nan),
        evg=np.where(g_ok, ev, np.inf),
        sg=np.where(g_ok, g_sig, nan))
    return maps, x, y, s_ok, g_ok
