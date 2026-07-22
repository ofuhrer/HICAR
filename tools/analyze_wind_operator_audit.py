#!/usr/bin/env python3
"""Analyze the first calibrated HICAR FGMRES Arnoldi relation.

The input is emitted by wind_iterative.F90 when
HICAR_WIND_OPERATOR_AUDIT is enabled.  The Arnoldi relation describes the
right-preconditioned operator A M^{-1}; it is therefore the relevant small
operator for deciding whether restarted/deflated GMRES can retain the slow
modes lost at an ordinary restart.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def load_audit(
    path: Path,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, int], dict[str, np.ndarray]]:
    metadata: dict[str, int] = {}
    h_entries: list[tuple[int, int, float]] = []
    projections: dict[str, dict[int, float]] = {"bootstrap": {}, "current": {}}
    spatial_basis_entries: list[tuple[int, int, float]] = []
    spatial_cell_counts: dict[int, float] = {}

    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            record = row["record"]
            name = row["name"]
            if record == "metadata":
                metadata[name] = int(row["value"])
            elif record == "hessenberg":
                h_entries.append((int(row["index_i"]), int(row["index_j"]), float(row["value"])))
            elif record == "projection":
                projections.setdefault(name, {})[int(row["index_i"])] = float(row["value"])
            elif record == "spatial" and name == "basis_sum":
                spatial_basis_entries.append(
                    (int(row["index_i"]), int(row["index_j"]), float(row["value"]))
                )
            elif record == "spatial" and name == "cell_count":
                spatial_cell_counts[int(row["index_j"])] = float(row["value"])

    dimension = metadata["arnoldi_dimension"]
    hbar = np.zeros((dimension + 1, dimension), dtype=np.float64)
    for i, j, value in h_entries:
        if i <= dimension + 1 and j <= dimension:
            hbar[i - 1, j - 1] = value

    dense_projections: dict[str, np.ndarray] = {}
    for name, values in projections.items():
        dense_projections[name] = np.array(
            [values.get(i, 0.0) for i in range(1, dimension + 1)], dtype=np.float64
        )
    spatial: dict[str, np.ndarray] = {}
    n_spatial_bins = (
        metadata.get("spatial_x_bins", 0)
        * metadata.get("spatial_y_bins", 0)
        * metadata.get("spatial_z_bins", 0)
    )
    if n_spatial_bins and spatial_basis_entries:
        basis_sums = np.zeros((dimension, n_spatial_bins), dtype=np.float64)
        for basis, bin_index, value in spatial_basis_entries:
            if basis <= dimension and bin_index <= n_spatial_bins:
                basis_sums[basis - 1, bin_index - 1] = value
        spatial["basis_sums"] = basis_sums
        spatial["cell_counts"] = np.array(
            [spatial_cell_counts.get(index, 0.0) for index in range(1, n_spatial_bins + 1)],
            dtype=np.float64,
        )
    return hbar, dense_projections, metadata, spatial


def complex_value(value: complex) -> dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag), "magnitude": float(abs(value))}


def normalized_projection(vector: np.ndarray, coordinates: np.ndarray) -> float:
    denominator = np.linalg.norm(vector) * np.linalg.norm(coordinates)
    if denominator == 0.0:
        return 0.0
    return float(abs(np.vdot(vector, coordinates)) / denominator)


def spatial_mode_summary(
    coefficients: np.ndarray,
    spatial: dict[str, np.ndarray],
    metadata: dict[str, int],
) -> dict | None:
    if "basis_sums" not in spatial:
        return None
    nx = metadata["spatial_x_bins"]
    ny = metadata["spatial_y_bins"]
    nz = metadata["spatial_z_bins"]
    counts = spatial["cell_counts"]
    valid = counts > 0.0
    coarse_sums = coefficients @ spatial["basis_sums"]
    coarse_means = np.zeros(coarse_sums.shape, dtype=np.complex128)
    coarse_means[valid] = coarse_sums[valid] / counts[valid]
    energy = counts * np.abs(coarse_means) ** 2
    captured = float(np.sum(energy))
    weights = energy / max(captured, np.finfo(float).tiny)

    indices = np.arange(counts.size)
    bx = indices % nx
    by = (indices // nx) % ny
    bz = indices // (nx * ny)
    centers = (
        (bx.astype(float) + 0.5) / nx,
        (by.astype(float) + 0.5) / ny,
        (bz.astype(float) + 0.5) / nz,
    )
    centroid = [float(np.sum(weights * coordinate)) for coordinate in centers]
    spread = [
        float(np.sqrt(np.sum(weights * (coordinate - center) ** 2)))
        for coordinate, center in zip(centers, centroid)
    ]
    dominant = np.argsort(energy)[::-1][: min(10, energy.size)]
    return {
        "grid_shape_xyz": [nx, ny, nz],
        "coarse_mean_energy_fraction": captured,
        "normalized_centroid_xyz": centroid,
        "normalized_spread_xyz": spread,
        "dominant_bins": [
            {
                "index_xyz": [int(bx[index]), int(by[index]), int(bz[index])],
                "fraction_of_coarse_energy": float(weights[index]),
                "mean_value": complex_value(complex(coarse_means[index])),
            }
            for index in dominant
            if valid[index]
        ],
    }


def analyze(
    hbar: np.ndarray,
    projections: dict[str, np.ndarray],
    metadata: dict[str, int],
    keep: int,
    spatial: dict[str, np.ndarray] | None = None,
) -> dict:
    dimension = hbar.shape[1]
    hessenberg = hbar[:dimension, :]
    h_last = float(hbar[dimension, dimension - 1])
    norm_h = float(np.linalg.norm(hessenberg, ord="fro"))

    commutator = hessenberg.T @ hessenberg - hessenberg @ hessenberg.T
    projected_nonnormality = float(
        np.linalg.norm(commutator, ord="fro") / max(norm_h * norm_h, np.finfo(float).tiny)
    )
    ordinary_values = np.linalg.eigvals(hessenberg)
    henrici_squared = max(norm_h * norm_h - float(np.sum(np.abs(ordinary_values) ** 2)), 0.0)
    henrici_departure = float(np.sqrt(henrici_squared) / max(norm_h, np.finfo(float).tiny))

    symmetric_part = 0.5 * (hessenberg + hessenberg.T)
    field_of_values_bounds = np.linalg.eigvalsh(symmetric_part)

    e_last = np.zeros(dimension, dtype=np.float64)
    e_last[-1] = 1.0
    transpose_condition = float(np.linalg.cond(hessenberg.T))
    try:
        correction = np.linalg.solve(hessenberg.T, e_last)
        harmonic_matrix = hessenberg + (h_last * h_last) * np.outer(correction, e_last)
        harmonic_values, harmonic_vectors = np.linalg.eig(harmonic_matrix)
        harmonic_status = "ok"
    except np.linalg.LinAlgError:
        harmonic_values, harmonic_vectors = np.linalg.eig(hessenberg)
        harmonic_status = "singular_hessenberg_fallback_to_ordinary_ritz"

    ordering = np.argsort(np.abs(harmonic_values))
    modes = []
    for order_index in ordering[: min(keep, dimension)]:
        vector = harmonic_vectors[:, order_index]
        coefficient_order = np.argsort(np.abs(vector))[::-1][: min(10, dimension)]
        mode = {
                "harmonic_ritz_value": complex_value(complex(harmonic_values[order_index])),
                "bootstrap_rhs_projection": normalized_projection(vector, projections["bootstrap"]),
                "current_rhs_projection": normalized_projection(vector, projections["current"]),
                "largest_basis_coefficients": [
                    {
                        "basis_index": int(index + 1),
                        "coefficient": complex_value(complex(vector[index])),
                    }
                    for index in coefficient_order
                ],
            }
        spatial_summary = spatial_mode_summary(vector, spatial or {}, metadata)
        if spatial_summary is not None:
            mode["spatial_structure"] = spatial_summary
        modes.append(mode)

    current = projections["current"]
    current_tail_fraction = float(
        np.linalg.norm(current[1:]) / max(np.linalg.norm(current), np.finfo(float).tiny)
    )
    return {
        "operator": "right_preconditioned_A_M_inverse",
        "arnoldi_dimension": dimension,
        "configured_restart": metadata.get("restart"),
        "bootstrap_rhs_saved": bool(metadata.get("bootstrap_rhs_saved", 0)),
        "hessenberg_transpose_condition_number": transpose_condition,
        "projected_nonnormality_commutator": projected_nonnormality,
        "henrici_departure": henrici_departure,
        "projected_field_of_values": {
            "minimum_real_rayleigh": float(field_of_values_bounds[0]),
            "maximum_real_rayleigh": float(field_of_values_bounds[-1]),
        },
        "current_residual_arnoldi_tail_fraction": current_tail_fraction,
        "harmonic_ritz_status": harmonic_status,
        "slow_harmonic_modes": modes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("audit_csv", type=Path)
    parser.add_argument("--keep", type=int, default=20, help="number of smallest-magnitude modes")
    parser.add_argument("--output", type=Path, help="write JSON here instead of stdout")
    args = parser.parse_args()

    hbar, projections, metadata, spatial = load_audit(args.audit_csv)
    result = analyze(hbar, projections, metadata, args.keep, spatial)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
