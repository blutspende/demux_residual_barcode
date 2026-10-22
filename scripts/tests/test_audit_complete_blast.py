import json
import subprocess
import sys
import unittest
from tempfile import NamedTemporaryFile, TemporaryDirectory
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_complete_blast import (
    barcode_locus, complete_identity, element_category, historical_identity, read_hsps, taxonomy_key,
)

SCRIPT = Path(__file__).resolve().parents[1] / "audit_complete_blast.py"


def hsp_line(qseqid, qlen, sstart, send, accession="ACC.1"):
    return (f"{qseqid}\t{accession}\t562\tN/A\t100.000\t{qlen}\t{qlen}\t100\t0\t0\t1\t{qlen}"
            f"\t{sstart}\t{send}\t1e-05\t48.1\tEscherichia coli plasmid pX\n")


class CompleteIdentityTests(unittest.TestCase):
    def setUp(self):
        self.record = dict(qlen=24, pident=100, mismatch=0, gapopen=0,
                           qcovs=100, length=24, qstart=1, qend=24,
                           sstart=101, send=124)

    def test_complete_both_coordinate_orientations(self):
        for query_start, query_end in [(1, 24), (24, 1)]:
            for subject_start, subject_end in [(101, 124), (124, 101)]:
                with self.subTest(query_start=query_start, subject_start=subject_start):
                    record = dict(self.record, qstart=query_start, qend=query_end,
                                  sstart=subject_start, send=subject_end)
                    self.assertTrue(complete_identity(record))

    def test_subject_coverage_does_not_prove_complete_hsp(self):
        record = dict(self.record, qstart=2, length=23, send=123)
        self.assertTrue(historical_identity(record))
        self.assertFalse(complete_identity(record))

    def test_each_constraint_is_required(self):
        for change in [dict(pident=99), dict(mismatch=1), dict(gapopen=1),
                       dict(gaps=1), dict(qlen=0), dict(qstart=2), dict(qend=23),
                       dict(length=25), dict(send=125), dict(sstart=0, send=23)]:
            with self.subTest(change=change):
                self.assertFalse(complete_identity(dict(self.record, **change)))

    def test_construct_query_length(self):
        self.assertTrue(complete_identity(dict(self.record, qlen=68, length=68,
                                              qend=68, send=168)))
        self.assertFalse(complete_identity(dict(self.record, qlen=68)))

    def test_query_coverage_is_not_the_controlling_field(self):
        self.assertTrue(complete_identity(dict(self.record, qcovs=99)))

    def test_exported_reverse_identifier(self):
        with NamedTemporaryFile(mode="w+", suffix=".tsv") as stream:
            stream.write("NB01_reverse\tACC.1\t123\tN/A\t100\t24\t24\t100\t0\t0"
                         "\t1\t24\t124\t101\t0.01\t48.1\tRecord title\n")
            stream.flush()
            record = next(read_hsps(stream.name))
            self.assertEqual(record["orientation"], "revcomp")
            self.assertTrue(complete_identity(record))

    def test_taxonomy_does_not_parse_titles_or_split_ambiguous_sets(self):
        self.assertEqual(taxonomy_key("123"), "taxid:123")
        self.assertEqual(taxonomy_key("456;123;123"), "taxid_set:123;456")
        self.assertEqual(taxonomy_key("N/A"), "unresolved")

    def test_category_precedence_and_neutral_genome_label(self):
        self.assertEqual(element_category("Phage plasmid complete genome"), "Plasmid")
        self.assertEqual(element_category("Eukaryotic chromosome"), "Chromosome/genome")
        self.assertEqual(element_category("Unknown record"), "Other/uncultured")


class ConstructPanelTests(unittest.TestCase):
    def test_identifiers_must_match_panel(self):
        with NamedTemporaryFile(mode="w+", suffix=".tsv") as stream:
            stream.write(hsp_line("NB01_revcomp_K46", 46, 146, 101))
            stream.flush()
            record = next(read_hsps(stream.name, "K46"))
            self.assertEqual((record["barcode"], record["orientation"]), ("NB01", "revcomp"))
            with self.assertRaises(ValueError):
                next(read_hsps(stream.name, "K40"))
            with self.assertRaises(ValueError):
                next(read_hsps(stream.name, "B24"))

    def test_barcode_locus_is_strand_and_orientation_invariant(self):
        loci = {
            barcode_locus(dict(barcode="NB01", saccver="ACC.1", orientation="forward",
                               qstart=1, qend=46, sstart=101, send=146), "K46"),
            barcode_locus(dict(barcode="NB01", saccver="ACC.1", orientation="revcomp",
                               qstart=1, qend=46, sstart=146, send=101), "K46"),
            barcode_locus(dict(barcode="NB01", saccver="ACC.1", orientation="forward",
                               qstart=1, qend=68, sstart=79, send=146), "K68"),
            barcode_locus(dict(barcode="NB01", saccver="ACC.1", orientation="revcomp",
                               qstart=1, qend=24, sstart=138, send=115), "B24"),
        }
        self.assertEqual(loci, {("NB01", "ACC.1", 115, 138)})

    def run_audit(self, directory, panels):
        arguments = [sys.executable, str(SCRIPT), "--out", str(Path(directory) / "out")]
        for panel, lines in panels.items():
            path = Path(directory) / f"{panel}.tsv"
            path.write_text("".join(lines))
            arguments += ["--panel", f"{panel}={path}"]
        result = subprocess.run(arguments, capture_output=True, text=True)
        report = json.loads((Path(directory) / "out" / "audit_summary.json").read_text())
        return result.returncode, report

    def test_nested_constructs_and_context(self):
        panels = {"B24": [hsp_line("NB01_forward", 24, 115, 138), hsp_line("NB01_reverse", 24, 138, 115),
                          hsp_line("NB01_forward", 24, 515, 538), hsp_line("NB01_reverse", 24, 538, 515)],
                  "K40": [hsp_line("NB01_forward_K40", 40, 107, 146),
                          hsp_line("NB01_revcomp_K40", 40, 146, 107)]}
        with TemporaryDirectory() as directory:
            code, report = self.run_audit(directory, panels)
            self.assertEqual(code, 0)
            checks = report["consistency"]
            self.assertTrue(checks["passed"])
            self.assertEqual(checks["nesting"]["K40_within_B24"], dict(outer_loci=1, missing_from_inner=0))
            self.assertEqual(checks["context"]["K40"]["loci"], 1)
            self.assertEqual(checks["context"]["B24"]["loci"], 2)
            context = (Path(directory) / "out" / "construct_context_pairs.csv").read_text().splitlines()
            self.assertEqual(context[1].split(",")[-3:], ["2", "1", "K40"])

    def test_construct_without_barcode_locus_fails(self):
        panels = {"B24": [hsp_line("NB01_forward", 24, 515, 538), hsp_line("NB01_reverse", 24, 538, 515)],
                  "K40": [hsp_line("NB01_forward_K40", 40, 107, 146),
                          hsp_line("NB01_revcomp_K40", 40, 146, 107)]}
        with TemporaryDirectory() as directory:
            code, report = self.run_audit(directory, panels)
            self.assertEqual(code, 1)
            self.assertEqual(report["consistency"]["nesting"]["K40_within_B24"]["missing_from_inner"], 1)

    def test_asymmetric_orientations_fail(self):
        panels = {"B24": [hsp_line("NB01_forward", 24, 115, 138)]}
        with TemporaryDirectory() as directory:
            code, report = self.run_audit(directory, panels)
            self.assertEqual(code, 1)
            self.assertFalse(report["consistency"]["symmetry"]["B24"]["orientation_locus_sets_equal"])


if __name__ == "__main__":
    unittest.main()