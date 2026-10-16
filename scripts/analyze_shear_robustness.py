"""Reanalyse the committed KiDS-1000 xi+ browser derivative; no network calls."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from research_analysis.shear_diagnostics import (  # noqa: E402
    covariance_audit, leave_one_angle_out, power_law_profile)


def analyze(payload: dict) -> dict:
    if payload.get('statistic') != 'xi_plus':
        raise ValueError('expected KiDS xi_plus data')
    results = []
    for pair in payload['pairs']:
        audit = covariance_audit(pair)
        leave_one_out = leave_one_angle_out(pair)
        profile = power_law_profile(pair)
        results.append({'pair_id': pair['id'],
                        'covariance': audit,
                        'leave_one_angle_out': leave_one_out,
                        'power_law': {
                            'best': profile['best'],
                            'pivot_arcmin': profile['pivot_arcmin'],
                            'nominal_residual_dof': profile['nominal_residual_dof'],
                            'profile': profile['profile'],
                            'warning': profile['warning']}})
    return {
        'source_fits_sha256': payload['source_fits_sha256'],
        'dataset': payload['dataset'],
        'n_pairs': len(results), 'pairs': results,
        'limitations': [
            'within-pair 9x9 covariance only; no cross-pair likelihood',
            'power-law is descriptive and does not infer cosmology',
            'leave-one-angle-out deltas are conditional and not independent',
            'no source redshift, intrinsic-alignment or scale-cut marginalization',
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=ROOT/'data/kids1000-xip.json')
    parser.add_argument('--output', type=Path, default=ROOT/'results/kids_covariance_robustness.json')
    args = parser.parse_args()
    output = analyze(json.loads(args.input.read_text(encoding='utf-8')))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(f"wrote {args.output} with {output['n_pairs']} within-pair diagnostics")

if __name__ == '__main__': main()
