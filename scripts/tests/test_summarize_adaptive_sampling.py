import csv
import importlib.util
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "summarize_adaptive_sampling.py"
SPEC = importlib.util.spec_from_file_location("summarize_adaptive_sampling", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

HEADER = ["read_id", "passes_filtering", "sequence_length_template", "end_reason", "mean_qscore_template",
          "barcode_arrangement"]
ROWS = [
    ("a", "TRUE", 400, "data_service_unblock_mux_change", 10.9),
    ("b", "TRUE", 500, "data_service_unblock_mux_change", 12.2),
    ("c", "TRUE", 900, "signal_positive", 14.5), ("d", "TRUE", 1200, "unblock_mux_change", 20.0),
    ("e", "TRUE", 1500, "signal_positive", 30.7), ("f", "FALSE", 50, "data_service_unblock_mux_change", 3.0),
]


class SummarizeAdaptiveSampling(unittest.TestCase):
    def test_median_of_counts(self):
        self.assertEqual(MODULE.median({5: 1}), 5)
        self.assertEqual(MODULE.median({1: 1, 3: 1}), 2)
        self.assertEqual(MODULE.median({1: 2, 9: 1}), 1)
        self.assertIsNone(MODULE.median({}))

    def test_nanostat_mean(self):
        self.assertAlmostEqual(MODULE.nanostat_mean({10: 1}), 10.0)
        self.assertAlmostEqual(MODULE.nanostat_mean({10: 1, 20: 1}), -10 * math.log10((0.1 + 0.01) / 2))
        self.assertIsNone(MODULE.nanostat_mean({}))

    def test_pass_reads_by_end_reason(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = Path(tmp) / "sequencing_summary.txt"
            summary.write_text("\t".join(HEADER) + "\n" +
                               "".join("\t".join(map(str, row)) + "\tbarcode01\n" for row in ROWS))
            out = Path(tmp) / "table.csv"
            subprocess.run([sys.executable, str(SCRIPT), "--summary", f"Exp 1={summary}", "--out", str(out)],
                           check=True)
            row = next(csv.DictReader(open(out)))
        # Pass-read Q-scores 10.9, 12.2, 14.5, 20.0 and 30.7; NanoStat truncates them to 10, 12, 14, 20 and 30.
        expected = -10 * math.log10(sum(10 ** (-q / 10) for q in (10, 12, 14, 20, 30)) / 5)
        self.assertEqual(row, {"experiment": "Exp 1", "pass_reads": "5", "rejected_pass_reads": "2",
                               "rejected_percent": "40.0", "median_length_rejected": "450.0",
                               "median_length_other_pass": "1200.0", "median_qscore_pass": "14.5",
                               "mean_qscore_pass": "17.66", "nanostat_mean_qscore_pass": str(round(expected, 4))})


if __name__ == "__main__":
    unittest.main()
