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


def load_audit(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, int]]:
    metadata: dict[str, int] = {}
    h_entries: list[tuple[int, int, float]] = []
    projections: dict[str, dict[int, float]] = {"bootstrap": {}, "current": {}}

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
    return hbar, dense_projections, metadata


def complex_value(value: complex) -> dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag), "magnitude": float(abs(value))}


def normalized_projection(vector: np.ndarray, coordinates: np.ndarray) -> float:
    denominator = np.linalg.norm(vector) * np.linalg.norm(coordinates)
    if denominator == 0.0:
        return 0.0
    return float(abs(np.vdot(vector, coordinates)) / denominator)


def analyze(hbar: np.ndarray, projections: dict[str, np.ndarray], metadata: dict[str, int], keep: int) -> dict:
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
        modes.append(
            {
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
        )

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

    hbar, projections, metadata = load_audit(args.audit_csv)
    result = analyze(hbar, projections, metadata, args.keep)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
