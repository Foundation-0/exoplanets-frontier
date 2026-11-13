#!/usr/bin/env python3
"""Two fixed synthetic controls of the original seeded BLS statistic, without detrending."""
from __future__ import annotations
import argparse
import ast
import hashlib
from importlib import metadata
import json
import math
from pathlib import Path
import platform
import numpy as np

SOURCE_HEAD = "883e28a7cecde05fdb6679bf574dbe78668332e7"
EXPECTED_SOURCE_HASHES = {
    "hz_core.py": "7977d141a473281c328c944cd5772670c205a059032e16af3e37a64e2e6c9a7e",
    "bls_search.py": "a5bb48b81db07d1d633b39c72cfea340e82e4f8c6f46b8ba8a755f64f3f40977",
}
SETTINGS = {
    'periods_days': [17.0, 25.0],
    'baseline_days': 80.0,
    'cadence_minutes': 30.0,
    'stellar_radius_solar': 0.55,
    'epoch_day': 5.0,
    'box_depth_ppm': 1000.0,
    'gaussian_noise_sigma_ppm': 100.0,
    'rng_seed': 314159,
    'duration_rule': 'unchanged bls_search.exp_dur(period, stellar_radius)',
    'noise_pairing': 'same fixed Gaussian noise vector for both periods',
    'sde_gate': 7.5,
}

def extract_function(root, file_name, function_names, env):
    tree = ast.parse((root / file_name).read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in function_names]
    assert {n.name for n in nodes} == set(function_names)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), file_name, 'exec'), env)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True, help='Path to the exoplanets-frontier checkout')
    root = parser.parse_args().repo.resolve()
    actual_hashes = {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                     for name in EXPECTED_SOURCE_HASHES}
    if actual_hashes != EXPECTED_SOURCE_HASHES:
        raise SystemExit('Source hashes differ from the reviewed revision. Control not run.')
    env = {'np': np, 'math': math}
    extract_function(root, 'bls_search.py', ['exp_dur', 'pont_beta'], env)
    extract_function(root, 'hz_core.py', ['coherent_seed_bls'], env)
    t = np.arange(0.0, SETTINGS['baseline_days'], SETTINGS['cadence_minutes'] / 1440.0)
    noise = np.random.default_rng(SETTINGS['rng_seed']).normal(
        0.0, SETTINGS['gaussian_noise_sigma_ppm'] * 1e-6, len(t))
    fe = np.full_like(t, SETTINGS['gaussian_noise_sigma_ppm'] * 1e-6)
    results = []
    for period in SETTINGS['periods_days']:
        duration = env['exp_dur'](period, SETTINGS['stellar_radius_solar'])
        phase = ((t - SETTINGS['epoch_day'] + 0.5 * period) % period) - 0.5 * period
        injected = np.abs(phase) <= 0.5 * duration
        f = 1.0 + noise - SETTINGS['box_depth_ppm'] * 1e-6 * injected
        pk = env['coherent_seed_bls'](t, f, fe, period, SETTINGS['stellar_radius_solar'])
        assert pk is not None, 'No seeded BLS result'
        results.append({
            'injected_period_days': period,
            'injected_duration_hours': duration * 24.0,
            'injected_in_transit_cadences': int(np.sum(injected)),
            'total_cadences': len(t),
            'fit_period_days': pk['p'],
            'relative_period_error': abs(pk['p'] - period) / period,
            'fit_epoch_day': pk['t0'],
            'fit_duration_hours': pk['dur'] * 24.0,
            'fit_depth_ppm': pk['dep'] * 1e6,
            'mes_raw': pk['mes'],
            'mes_corrected': pk['mes_c'],
            'sde': pk['sde'],
            'n_transits': pk['ntr'],
            'would_fail_sde_gate': pk['sde'] < SETTINGS['sde_gate'],
        })
    report = {
        'settings': SETTINGS,
        'python_version': platform.python_version(),
        'package_versions': {name: metadata.version(name) for name in ['numpy', 'astropy', 'pyerfa', 'astropy-iers-data', 'PyYAML', 'packaging']},
        'reviewed_source_head': SOURCE_HEAD,
        'source_hashes_match_reviewed_revision': True,
        'source_files_sha256': actual_hashes,
        'control_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Seeded BLS statistic only. No detrending, archive data, full vetter, stellar variability, completeness estimate, or target disposition.',
        'results': results,
    }
    print(json.dumps(report, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
