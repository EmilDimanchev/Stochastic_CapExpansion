"""
Independent check of a cut archive, read the way the dashboard's data step would read it.

Usage:
  python3 tests/check_cut_archive_h5py.py outputs/Test_cx_1/Summary/cut_archive_pareto.h5 \
      outputs/Test_cx_1/Summary/Summary_pareto.csv

It shares no code with src/results/Export_cuts.jl, so a wrong column order, scenario order or
array orientation in the file makes it fail. For every plan and scenario, the floor from the
file, LB_s(x) = max_k (alpha_k + beta_k . x), is compared with the operating cost the run wrote
to the Summary CSV.
"""

import sys

import h5py
import numpy as np
import pandas as pd


def text(values):
    """h5py returns strings as bytes; decode them."""
    return [v.decode() if isinstance(v, bytes) else v for v in values]


def read_beta(f, n_vars):
    """beta as a dense K x n array, from either storage format."""
    if f["meta/beta_format"][()].decode() == "csr":
        group = f["cuts/beta_csr"]
        indptr, indices, data = group["indptr"][()], group["indices"][()], group["data"][()]
        beta = np.zeros((len(indptr) - 1, n_vars))
        for k in range(len(indptr) - 1):
            entries = range(indptr[k], indptr[k + 1])
            beta[k, indices[entries.start : entries.stop]] = data[entries.start : entries.stop]
        return beta
    return f["cuts/beta"][()]


def main(archive_path, summary_path, relative_tolerance=1e-4):
    summary = pd.read_csv(summary_path).set_index("Index")

    with h5py.File(archive_path, "r") as f:
        names = text(f["vars/name"][()])
        plan_ids = text(f["plans/id"][()])
        plan_x = f["plans/x"][()].astype(np.float64)  # N x n
        scenario = f["cuts/scenario"][()]  # 1-based
        alpha = f["cuts/alpha"][()]
        beta = read_beta(f, len(names)).astype(np.float64)  # K x n
        labels = text(f["meta/scenario_label"][()])
        probabilities = f["meta/scenario_probability"][()]
        failed = set(text(f["meta/failed_plans"][()]))

    assert beta.shape == (len(alpha), len(names)), f"beta has shape {beta.shape}"
    assert plan_x.shape == (len(plan_ids), len(names)), f"plans/x has shape {plan_x.shape}"
    assert abs(probabilities.sum() - 1) < 1e-9, "scenario probabilities don't sum to 1"

    # plans/x matches the Summary's capacity columns, column by column
    summary_x = summary.loc[plan_ids, names].to_numpy()
    x_error = np.abs(plan_x - summary_x) / np.maximum(np.abs(summary_x), 1.0)
    print(f"plans/x vs Summary capacity columns: max relative difference {x_error.max():.2e}")
    assert x_error.max() < 1e-5, "plans/x doesn't match the Summary"

    # Floor per plan and scenario, against the Summary's per-scenario operating cost
    checked_plans = [i for i, plan in enumerate(plan_ids) if plan not in failed]
    gaps = []
    for s, label in enumerate(labels, start=1):
        rows = scenario == s
        assert rows.any(), f"scenario {s} has no cuts"
        floor = (alpha[rows, None] + beta[rows] @ plan_x[checked_plans].T).max(axis=0)
        cost = summary.loc[[plan_ids[i] for i in checked_plans], f"OperatingCost_{label}"].to_numpy()
        gaps.append((cost - floor) / np.maximum(np.abs(cost), 1.0))
    gaps = np.concatenate(gaps)

    print(f"Cuts {len(alpha)}, scenarios {len(labels)}, variables {len(names)}, plans {len(plan_ids)}")
    print(f"(Op_s - LB_s) / Op_s over {gaps.size} pairs: min {gaps.min():.2e}, max {gaps.max():.2e}, mean {gaps.mean():.2e}")
    assert np.abs(gaps).max() < relative_tolerance, "floors from the file don't match the Summary's op costs"
    print("PASSED")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
