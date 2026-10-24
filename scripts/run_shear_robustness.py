"""Run KiDS-1000 within-pair checks on the repository's committed JSON."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from shear_robustness import covariance_diagnostic, leave_one_angle_out, null_diagnostic, powerlaw_profile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parents[1] / "data/kids1000-xip.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.data.read_text(encoding="utf-8"))
    if payload.get("statistic") != "xi_plus" or len(payload.get("pairs", [])) != 15:
        raise ValueError("expected 15 KiDS xi+ tomographic pairs")
    results = []
    for pair in payload["pairs"]:
        results.append({"pair": pair["id"], "null": null_diagnostic(pair),
                        "covariance": covariance_diagnostic(pair),
                        "leave_one_angle_out": leave_one_angle_out(pair),
                        "empirical_powerlaw": powerlaw_profile(pair)})
    report = {"dataset": payload["dataset"], "source_fits_sha256": payload["source_fits_sha256"],
              "limitations": "Within-pair covariance only; no inter-pair joint significance or cosmological inference.",
              "pairs": results}
    content = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(content, encoding="utf-8")
    else:
        print(content)


if __name__ == "__main__":
    main()
