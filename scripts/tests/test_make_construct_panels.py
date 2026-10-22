import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from make_construct_panels import (
    PANELS, SSP_22, barcode_coordinates, build, construct, panel_queries, read_barcodes,
    reverse_complement,
)

NB01 = "CACAAAGACACCGACAACTTTCTT"
NB03 = "CCTGGTAACTGGGACACAAGACTC"
PROJECT_BARCODES = Path(__file__).resolve().parents[2] / "barcodes.csv"
EXPECTED_COORDINATES = {"B24": ((1, 24), (1, 24)), "K40": ((9, 32), (9, 32)),
                        "K46": ((15, 38), (9, 32)), "K68": ((37, 60), (9, 32))}


def write_csv(path, rows):
    with open(path, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["Component", "Forward sequence", "Reverse sequence"])
        writer.writerows(rows)


class ConstructTests(unittest.TestCase):
    def test_nb01_constructs(self):
        self.assertEqual(construct("B24", NB01), NB01)
        self.assertEqual(construct("K40", NB01), "AAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT")
        self.assertEqual(construct("K46", NB01), "ATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT")
        self.assertEqual(construct("K68", NB01),
                         "CCTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT")

    def test_adapter_and_flank_join_once(self):
        sequence = construct("K68", NB01)
        self.assertEqual(sequence.count("ATTGC"), 1)
        self.assertNotIn("ATTGCATTGC", sequence)

    def test_lengths_and_barcode_coordinates(self):
        for panel, lengths in zip(PANELS, (24, 40, 46, 68)):
            queries = list(panel_queries(panel, [("NB01", NB01, reverse_complement(NB01))]))
            for (_, _, orientation, sequence), expected in zip(queries, EXPECTED_COORDINATES[panel]):
                with self.subTest(panel=panel, orientation=orientation):
                    self.assertEqual(len(sequence), lengths)
                    start, end = barcode_coordinates(panel, orientation)
                    self.assertEqual((start, end), expected)
                    barcode = NB01 if orientation == "forward" else reverse_complement(NB01)
                    self.assertEqual(sequence[start - 1:end], barcode)

    def test_second_query_is_whole_construct_reverse_complement(self):
        for panel in PANELS:
            forward, revcomp = panel_queries(panel, [("NB01", NB01, reverse_complement(NB01))])
            self.assertEqual(revcomp[3], reverse_complement(forward[3]))
        self.assertEqual([query[0] for query in panel_queries("K46", [("NB01", NB01, "")])],
                         ["NB01_forward_K46", "NB01_revcomp_K46"])

    def test_nesting_and_absence_of_cdna_primer(self):
        for barcode in (NB01, NB03):
            sequences = [construct(panel, barcode) for panel in PANELS]
            for inner, outer in zip(sequences, sequences[1:]):
                self.assertIn(inner, outer)
            for sequence in sequences:
                self.assertNotIn(SSP_22, sequence)
                self.assertNotIn(reverse_complement(SSP_22), sequence)


class InputTests(unittest.TestCase):
    def test_rejects_reverse_column_that_is_not_reverse_complement(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "barcodes.csv"
            write_csv(path, [["NB01", NB01, NB01]])
            with self.assertRaises(SystemExit):
                read_barcodes(path)

    def test_subset_and_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "barcodes.csv"
            write_csv(path, [["NB01", NB01, reverse_complement(NB01)],
                             ["NB03", NB03, reverse_complement(NB03)]])
            manifest = build(path, Path(directory) / "queries", subset={"NB03"})
            self.assertEqual(manifest["barcodes"], ["NB03"])
            fasta = (Path(directory) / "queries" / "K68.fasta").read_text()
            self.assertEqual(fasta.count(">"), 2)
            self.assertIn(">NB03_revcomp_K68\n" + reverse_complement(construct("K68", NB03)) + "\n", fasta)
            with open(Path(directory) / "queries" / "query_map.tsv") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(rows), 8)
            saved = json.loads((Path(directory) / "queries" / "panel_manifest.json").read_text())
            self.assertEqual(saved["panels"]["K46"]["barcode_coordinates"]["forward"], [15, 38])


@unittest.skipUnless(PROJECT_BARCODES.exists(), "project barcodes.csv not available")
class ProjectBarcodeTests(unittest.TestCase):
    def test_all_96_barcodes_and_legacy_barcode_only_format(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = build(PROJECT_BARCODES, Path(directory))
            self.assertEqual(len(manifest["barcodes"]), 96)
            for panel in PANELS:
                text = (Path(directory) / f"{panel}.fasta").read_text()
                identifiers = [line for line in text.splitlines() if line.startswith(">")]
                self.assertEqual(len(identifiers), 192)
                self.assertEqual(len(set(identifiers)), 192)
            legacy = []
            with open(PROJECT_BARCODES, newline="") as stream:
                for row in csv.DictReader(stream):
                    legacy.append(f">{row['Component'].strip()}_forward\n{row['Forward sequence'].strip()}\n")
                    legacy.append(f">{row['Component'].strip()}_reverse\n{row['Reverse sequence'].strip()}\n")
            self.assertEqual((Path(directory) / "B24.fasta").read_text(), "".join(legacy))


if __name__ == "__main__":
    unittest.main()
