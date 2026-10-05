"""Derive the compact browser dataset from the official KiDS-1000 2PCF FITS file.

The output preserves every xi+ tomographic pair (15 pairs x 9 angular bins),
diagonal uncertainties, the full within-pair covariance, and the five source
redshift distributions. Detection statistics test the zero-shear null with
the corresponding 9x9 covariance submatrix; they are not cosmological fits.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "research" / "source" / "kids_xipm.fits"
OUTPUT = ROOT / "data" / "kids1000-xip.json"
SOURCE_URL = "https://kids.strw.leidenuniv.nl/DR4/data_files/KiDS1000_cosmic_shear_data_release.tgz"
ARCHIVE_SHA256 = "97729199aeb23238767921f58f36341317aa5251beab99e2575eb35fe78c7bb3"


def numbers(values):
    return [float(value) for value in values]


def main():
    source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    with fits.open(SOURCE) as hdus:
        table = hdus["xiP"].data
        covariance = np.asarray(hdus["COVMAT"].data[: len(table), : len(table)], dtype=float)
        nz = hdus["NZ_SOURCE"].data
        pairs = []
        for bin1 in range(1, 6):
            for bin2 in range(bin1, 6):
                indices = np.where((table["BIN1"] == bin1) & (table["BIN2"] == bin2))[0]
                vector = np.asarray(table["VALUE"][indices], dtype=float)
                block = covariance[np.ix_(indices, indices)]
                chi_square = float(vector @ np.linalg.solve(block, vector))
                pairs.append({
                    "id": f"{bin1}x{bin2}",
                    "bin1": bin1,
                    "bin2": bin2,
                    "theta_arcmin": numbers(table["ANG"][indices]),
                    "xi_plus": numbers(vector),
                    "sigma": numbers(np.sqrt(np.diag(block))),
                    "covariance": [numbers(row) for row in block],
                    "null_chi2": chi_square,
                    "degrees_of_freedom": int(len(indices)),
                    "null_p_value": float(chi2.sf(chi_square, len(indices))),
                    "correlated_snr": float(np.sqrt(chi_square)),
                })

        z_mid = np.asarray(nz["Z_MID"], dtype=float)
        redshift_bins = []
        for index in range(1, 6):
            density = np.asarray(nz[f"BIN{index}"], dtype=float)
            norm = float(np.trapz(density, z_mid))
            redshift_bins.append({
                "bin": index,
                "mean_z": float(np.trapz(z_mid * density, z_mid) / norm),
                "z_mid": numbers(z_mid),
                "density": numbers(density / norm),
            })

    payload = {
        "dataset": "KiDS-1000 cosmic shear two-point correlation functions",
        "statistic": "xi_plus",
        "release_page": "https://kids.strw.leidenuniv.nl/DR4/KiDS-1000_cosmicshear.php",
        "source_archive_url": SOURCE_URL,
        "source_archive_sha256": ARCHIVE_SHA256,
        "source_fits_sha256": source_hash,
        "source_fits_hdu": "xiP",
        "citation": "Asgari et al. 2021, Astronomy & Astrophysics 645, A104",
        "scope_note": "Observed KiDS-1000 xi+ data. Null statistics test xi+=0 for one bin pair using its 9x9 covariance block; they are not cosmological parameter constraints.",
        "survey_area_deg2": 1006,
        "source_bins": 5,
        "angular_bins_per_pair": 9,
        "pairs": pairs,
        "redshift_distributions": redshift_bins,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT} with {len(pairs)} tomographic pairs")


if __name__ == "__main__":
    main()
