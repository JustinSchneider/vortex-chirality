"""MaNGA DR17 access for the Gate G1 MaNGA test (PREREG amendments c-e).

Parent sample from drpall/dapall, and per-galaxy DAP maps (HYB10-MILESHC-
MASTARSSP) fetched from the SDSS Marvin API with verified TLS and cached as
.npz files. Nothing here prints or summarises rotation amplitudes; selection
code consumes velocity maps only through position-angle fits.
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
MANGA = ROOT / "data" / "manga"
CACHE = MANGA / "maps"
SAS = ("https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/"
       "HYB10-MILESHC-MASTARSSP/{plate}/{ifu}/"
       "manga-{pifu}-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz")
HA_CHANNEL = 23            # 'Ha-6564' is C24 in the EMLINE_* extensions
SIGMACORR_CHANNEL = 0      # DR17 guidance: use channel C1 'resolution difference'
API = "https://magrathea.sdss.org/marvin/api/maps/{pifu}/HYB10/MILESHC-MASTARSSP/map/{prop}/{chan}/"
SPAXEL_ARCSEC = 0.5
Q0 = 0.2                      # intrinsic thickness for inclinations
MAIN_BITS = (1 << 10) | (1 << 11) | (1 << 12)   # primary, secondary, color-enhanced
CRITICAL_BIT = 1 << 30

# (property, channel) for each map used.
MAPS = {
    "stellar_vel": ("stellar_vel", "None"),
    "stellar_sigma": ("stellar_sigma", "None"),
    "ha_vel": ("emline_gvel", "ha_6564"),
    "ha_sigma": ("emline_gsigma", "ha_6564"),
    "stellar_sigmacorr": ("stellar_sigmacorr", "resolution_difference"),
    "ha_instsigma": ("emline_instsigma", "ha_6564"),
    "ha_flux": ("emline_gflux", "ha_6564"),
}


def inclination_deg(ba, q0=Q0):
    cos2 = (np.asarray(ba, float) ** 2 - q0 ** 2) / (1 - q0 ** 2)
    return np.degrees(np.arccos(np.sqrt(np.clip(cos2, 0, 1))))


def parent_sample() -> pd.DataFrame:
    """Main-sample galaxies with good DRP and DAP quality, one row per galaxy."""
    drp = fits.getdata(MANGA / "drpall-v3_1_1.fits", "MANGA")
    dap = fits.getdata(MANGA / "dapall-v3_1_1-3.1.0.fits",
                       "HYB10-MILESHC-MASTARSSP")
    d = pd.DataFrame({
        "plateifu": drp["plateifu"].astype(str),
        "mangaid": drp["mangaid"].astype(str),
        "mngtarg1": drp["mngtarg1"].astype(np.int64),
        "drp3qual": drp["drp3qual"].astype(np.int64),
        "z": drp["nsa_z"].astype(float),
        "logMstar": np.log10(np.clip(drp["nsa_elpetro_mass"].astype(float), 1, None)),
        "ba": drp["nsa_elpetro_ba"].astype(float),
        "phi": drp["nsa_elpetro_phi"].astype(float) if "nsa_elpetro_phi" in drp.columns.names else np.nan,
        "Re_arcsec": drp["nsa_elpetro_th50_r"].astype(float),
        "sersic_n": drp["nsa_sersic_n"].astype(float),
    })
    a = pd.DataFrame({
        "plateifu": dap["PLATEIFU"].astype(str),
        "dapdone": dap["DAPDONE"].astype(bool),
        "dapqual": dap["DAPQUAL"].astype(np.int64),
        "adist_mpc": dap["ADIST_Z"].astype(float),
    })
    d = d.merge(a, on="plateifu", how="inner")
    keep = ((d["mngtarg1"] & MAIN_BITS) != 0) \
        & ((d["drp3qual"] & CRITICAL_BIT) == 0) \
        & d["dapdone"] & ((d["dapqual"] & CRITICAL_BIT) == 0) \
        & (d["ba"] > 0) & (d["Re_arcsec"] > 0)
    d = d[keep].copy()
    # One observation per galaxy: keep the first (repeat observations rare).
    d = d.drop_duplicates("mangaid", keep="first")
    d["inc"] = inclination_deg(d["ba"])
    d["kpc_per_arcsec"] = d["adist_mpc"] * 1e3 * np.pi / (180 * 3600)
    return d.set_index("plateifu")


def _sas_extract(pifu, session):
    """Download one MAPS file, cache every map in MAPS as .npz, delete it."""
    import tempfile
    plate, ifu = pifu.split("-")
    url = SAS.format(plate=plate, ifu=ifu, pifu=pifu)
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "maps.fits.gz"
        for attempt in range(5):
            try:
                with session.get(url, stream=True, timeout=120) as r:
                    r.raise_for_status()
                    with open(f, "wb") as fh:
                        for chunk in r.iter_content(1 << 20):
                            fh.write(chunk)
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(10 * (attempt + 1))
        with fits.open(f) as h:
            def ext(name, ch=None):
                d = h[name].data
                return d if ch is None else d[ch]

            def triple(base, ch=None):
                return {"value": ext(base, ch), "ivar": ext(base + "_IVAR", ch),
                        "mask": ext(base + "_MASK", ch)}

            def single(name, ch):
                v = ext(name, ch)
                return {"value": v, "ivar": np.ones_like(v),
                        "mask": np.zeros(v.shape, dtype=np.int64)}

            maps = {
                "stellar_vel": triple("STELLAR_VEL"),
                "stellar_sigma": triple("STELLAR_SIGMA"),
                "stellar_sigmacorr": single("STELLAR_SIGMACORR", SIGMACORR_CHANNEL),
                "ha_vel": triple("EMLINE_GVEL", HA_CHANNEL),
                "ha_sigma": triple("EMLINE_GSIGMA", HA_CHANNEL),
                "ha_instsigma": single("EMLINE_INSTSIGMA", HA_CHANNEL),
                "ha_flux": triple("EMLINE_GFLUX", HA_CHANNEL),
            }
            for key, m in maps.items():
                np.savez_compressed(CACHE / f"{pifu}_{key}.npz",
                                    **{k: np.asarray(v, dtype=float)
                                       for k, v in m.items()})


def fetch_maps(pifu, keys=("stellar_vel", "ha_vel", "ha_flux"), session=None):
    """Return cached maps, filling the cache from SAS MAPS files if needed."""
    CACHE.mkdir(parents=True, exist_ok=True)
    s = session or requests.Session()
    if not all((CACHE / f"{pifu}_{k}.npz").exists() for k in keys):
        _sas_extract(pifu, s)
    out = {}
    for key in keys:
        z = np.load(CACHE / f"{pifu}_{key}.npz")
        out[key] = {k: z[k] for k in z.files}
    return out


def fetch_maps_marvin(pifu, keys=("stellar_vel", "ha_vel", "ha_flux"), pause=1.0,
                      session=None):
    """Marvin API route (rate-limited; kept for reference, not used)."""
    """Fetch maps (value, ivar, mask) via Marvin, caching to data/manga/maps."""
    CACHE.mkdir(parents=True, exist_ok=True)
    out = {}
    s = session or requests.Session()
    for key in keys:
        f = CACHE / f"{pifu}_{key}.npz"
        if f.exists():
            z = np.load(f)
            out[key] = {k: z[k] for k in z.files}
            continue
        prop, chan = MAPS[key]
        for attempt in range(8):
            try:
                r = s.post(API.format(pifu=pifu, prop=prop, chan=chan),
                           data={"release": "DR17"}, timeout=90)
                if r.status_code == 429:
                    # Rate limited: honour Retry-After, else back off.
                    wait = float(r.headers.get("Retry-After", 0) or 0)
                    time.sleep(max(wait, 15 * 2 ** min(attempt, 4)))
                    continue
                r.raise_for_status()
                dd = r.json()
                if dd.get("error"):
                    raise RuntimeError(dd["error"])
                m = {k: np.asarray(dd["data"][k], dtype=float)
                     for k in ("value", "ivar", "mask")}
                np.savez_compressed(f, **m)
                out[key] = m
                break
            except Exception:
                if attempt == 7:
                    raise
                time.sleep(5 * (attempt + 1))
        else:
            raise RuntimeError(f"rate limited: {pifu} {key}")
        time.sleep(pause)
    return out
