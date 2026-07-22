from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "analyze_wind_operator_audit.py"
SPEC = importlib.util.spec_from_file_location("analyze_wind_operator_audit", MODULE_PATH)
ANALYZER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(ANALYZER)


class SpatialAuditTests(unittest.TestCase):
    def test_spatial_mode_summary_reports_location_and_captured_energy(self) -> None:
        spatial = {
            "basis_sums": np.array([[2.0, 0.0], [0.0, 3.0]]),
            "cell_counts": np.array([2.0, 3.0]),
        }
        metadata = {"spatial_x_bins": 2, "spatial_y_bins": 1, "spatial_z_bins": 1}
        summary = ANALYZER.spatial_mode_summary(np.array([1.0, 0.0]), spatial, metadata)
        self.assertIsNotNone(summary)
        self.assertAlmostEqual(summary["coarse_mean_energy_fraction"], 2.0)
        self.assertEqual(summary["dominant_bins"][0]["index_xyz"], [0, 0, 0])
        self.assertAlmostEqual(summary["normalized_centroid_xyz"][0], 0.25)

    def test_loader_accepts_spatial_records(self) -> None:
        rows = [
            ("metadata", "arnoldi_dimension", "", "", 1),
            ("metadata", "restart", "", "", 1),
            ("metadata", "bootstrap_rhs_saved", "", "", 0),
            ("metadata", "spatial_x_bins", "", "", 1),
            ("metadata", "spatial_y_bins", "", "", 1),
            ("metadata", "spatial_z_bins", "", "", 1),
            ("hessenberg", "H", 1, 1, 2.0),
            ("hessenberg", "H", 2, 1, 0.1),
            ("projection", "current", 1, "", 1.0),
            ("spatial", "cell_count", 0, 1, 4.0),
            ("spatial", "basis_sum", 1, 1, 2.0),
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "audit.csv"
            with path.open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("record", "name", "index_i", "index_j", "value"))
                writer.writerows(rows)
            hbar, projections, metadata, spatial = ANALYZER.load_audit(path)
        self.assertEqual(hbar.shape, (2, 1))
        self.assertEqual(projections["current"][0], 1.0)
        self.assertEqual(metadata["spatial_x_bins"], 1)
        self.assertEqual(spatial["basis_sums"][0, 0], 2.0)


if __name__ == "__main__":
    unittest.main()
