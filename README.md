# Demultiplexing and Residual Barcode Analysis

Calculation scripts for comparing Dorado and Barbell preprocessing, auditing complete-query barcode-reference matches, and checking detector sensitivity. The main experimental workflow uses the ONT SQK-NBD114.96 Native Barcoding Kit (kit name `SQK-NBD114-96` in Dorado and Barbell); the reference screen is a focused case study, not a general database-contamination census or cross-platform benchmark.

## Privacy and Repository Contents

This repository contains calculation/pipeline scripts, regression tests and this README only. It contains no sequencing inputs, donor or sample records, read identifiers, per-experiment measurements, result tables, taxonomy responses, manuscript files or reviews. Users supply their own data and barcode definitions locally.

The pipeline can produce sensitive intermediate files, including FASTQ output, read-level SQLite databases, examples of problematic reads and per-experiment reports. Those local outputs are **not approved for publication** merely because the code is public. Keep them outside Git and follow the applicable data-privacy policy. Do not use `git add .` in a directory containing study data.

A separately distributed `supplementary_source_data.zip` can contain the same release scripts and selected highly processed, cross-experiment summaries. That archive is not stored in this repository. A Zenodo version DOI should identify an archived release of this code; no DOI is asserted here before it has been issued.

## Requirements

- Python 3.12 or later, with `regex`, `numpy`, `scipy`, `pandas` and `matplotlib`.
- Dorado 1.1.1 and Dorado 2.0.0 for the two basecalling configurations. The pipeline uses a `dorado200` conda environment for the latter by default; override with `--conda-env-200`.
- Barbell 0.3.2, samtools 1.19.2 and NanoStat 1.6.0 for the recorded preprocessing workflow. Supply executable paths when they are not on PATH.
- ONT POD5 tooling if FAST5 conversion is needed; the pipeline takes POD5 input.
- Optional BLAST+ 2.15.0 (`blastn`, `blastdbcmd`) and a locally prepared nucleotide database for barcode screening. The construct cross-check additionally needs network access to NCBI E-utilities.
- Optional R with `ggplot2`, `readr` and Cairo graphics support for the current GraphPad-style Figure 2 renderer. It generates PNG/PDF files, not GraphPad project files.
- Optional minimap2 only for the legacy auxiliary barcode scan; it is not the approximate-regex endpoint described below.

Use an isolated environment:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install regex numpy scipy pandas matplotlib
```

The supported CLI entry points provide `--help`. Test the external tools independently before a full run. Thread counts and GPU requirements should match the host; the recorded 120-thread setting is not a recommended minimum.

## Preprocessing Conditions

| Condition | Configuration |
| --- | --- |
| C0 | Dorado 1.1.1, HAC v5.0.0, no demultiplexing or trimming |
| C1 | Dorado 2.0.0, HAC v6.0.0, no demultiplexing or trimming |
| C2/C3 | Corresponding Dorado demultiplexing and trimming |
| C4/C5 | Barbell kit mode with its default filter pattern set, applied to C0/C1 |
| C6 | Barbell (default pattern set) applied to the concatenated C3 output; a selected diagnostic subset |

The release and basecalling model change together. C4/C5 are parallel alternatives to C2/C3, not sequential processing of them. C6 is not a matched whole-population comparator.

The reported results use `barbell_paper.py --barbell-mode default`, which runs `barbell kit` without `--maximize`; the Barbell logs then report the kit type `Safe`. In Barbell 0.3.2, annotation does not depend on the pattern set, which acts only in the filter step. For this kit, the default patterns accept a barcode at either read end alone or concordant barcodes at both ends. The script default, `--barbell-mode maximize`, adds `--maximize` and was not used for the reported results. Retention and residual-pattern removal do not establish correct sample assignment or clinical safety.

## Pipeline

Run from the repository root. The barcode CSV is supplied by the user and must have columns `Component`, `Forward sequence` and `Reverse sequence`. No barcode-input file is bundled.

```bash
.venv/bin/python scripts/barbell_paper.py \
  --pod5-dir /path/to/retained/pod5 \
  --out outputs/experiment_a_demux \
  --barcode-csv /path/to/barcodes.csv \
  --dorado-111 /path/to/dorado-1.1.1/bin/dorado \
  --conda-env-200 dorado200 \
  --barbell /path/to/barbell \
  --barbell-mode default \
  --samtools /path/to/samtools \
  --python .venv/bin/python \
  --threads 8 \
  --skip-compare --skip-figures
```

Paths in this release are portable defaults or explicit arguments; they do not encode the source study's experiment labels. Input selection, including adaptive-sampling pass/fail handling, is performed before invoking the pipeline and must be appropriate to the user's data. With `--barbell-mode default`, the C4-C6 outputs are written to `exp_c4_barbell111_default`, `exp_c5_barbell200_default` and `exp_c6_dorado200_barbell_default`.

`barbell_paper.py` does not pass `--max-edits` to the scanner, so the residual-barcode scan of the reported results was run separately on all seven conditions:

```bash
D=outputs/experiment_a_demux
.venv/bin/python scripts/compare_dorado_barbell_outputs.py \
  --condition-dirs $D/exp_c0_dorado111_raw $D/exp_c1_dorado200_raw \
                   $D/exp_c2_dorado111_demux $D/exp_c3_dorado200_demux \
                   $D/exp_c4_barbell111_default/barbell_demux \
                   $D/exp_c5_barbell200_default/barbell_demux \
                   $D/exp_c6_dorado200_barbell_default/barbell_demux \
  --condition-names c0 c1 c2 c3 c4 c5 c6 \
  --barcode-csv /path/to/barcodes.csv \
  --out $D/cmp_all \
  --detector occurrences --max-edits 2 \
  --decoys 2 --decoy-seed 20260929 \
  --skip-nanostat --skip-barcode-scan \
  --threads 8
```

The study used `--decoys 2` for three of its eight experiments and `--decoys 0` for the others (see [Residual Barcode Detector](#residual-barcode-detector)). `--skip-barcode-scan` skips the legacy **minimap2** auxiliary scan, not the regex endpoint, and `--skip-nanostat` skips NanoStat QC; omit it and pass `--nanostat` to write per-condition NanoStat summaries. The regex scan is enabled by `--barcode-csv`; omitting that CSV skips residual-pattern detection. Choose a new output directory for a new analysis. `--overwrite` replaces an existing comparison database and should be used only deliberately.

## Aggregation and Paired Statistics

The local directory pattern is `outputs/{experiment}_demux/cmp_all/`. Supply the experiment names explicitly; the examples below are placeholders, not subject identifiers.

```bash
.venv/bin/python scripts/barbell_figures_summa.py \
  --base-dir outputs --experiments experiment_a experiment_b \
  --out figures/summary_figures --skip-figures

.venv/bin/python scripts/barbell_figure_f_paired.py \
  --base-dir outputs --experiments experiment_a experiment_b \
  --out figures/summary_figures/figure_f_violin

.venv/bin/python scripts/make_barcode_multiplicity_csvs.py \
  --base-dir outputs --experiments experiment_a experiment_b \
  --out figures/summary_figures/graphpad
```

`--base-dir` selects the result tree; `barbell_figure_f_paired.py` and `make_barcode_multiplicity_csvs.py` also take `--cmp-subdir` for a comparison folder other than `cmp_all`. `make_graphpad_fig1e_table.py` writes the per-experiment barcode-positive table from `summary_e_barcode_pct.csv`.

`barbell_figures.py` generates per-condition summaries from one comparison database. `barbell_figure_f_violin.py` calculates location-rate summaries and optional rarefaction displays from supplied experiment directories. `make_summary_f_log.py` takes `--input`, `--experiments` and `--out` to generate log-scale location summaries. Their generated replicate-level files are local analysis outputs, not privacy-filtered publication data.

The primary inferential endpoint is the paired, experiment-level internal-match rate in reads per million output reads. The paired script reports unadjusted Wilcoxon signed-rank tests and ratios. Its legacy plotting defaults are not a claim about the formatting of any final manuscript figure.

## Barcode-Database Screen

The BLAST driver builds four query panels from the barcode CSV with `make_construct_panels.py` (sequences from Dorado v1.1.1) and searches each against a user-supplied database prefix. Each panel has two queries per barcode: the construct and its whole-construct reverse complement.

| Panel | Query | Length |
|---|---|---|
| B24 | barcode | 24 |
| K40 | `AAGGTTAA` + barcode + `CAGCACCT` (native ligation flanks) | 40 |
| K46 | `ATTGCTAAGGTTAA` + barcode + `CAGCACCT` (Dorado `kit_nb96` search window) | 46 |
| K68 | `CCTGTACTTCGTTCAGTTACGTATTGC` (LSK110 adapter front) + `TAAGGTTAA` + barcode + `CAGCACCT` | 68 |

```bash
BLAST_BIN=/path/to/ncbi-blast-2.15.0+/bin THREADS=8 \
  bash scripts/run_barcode_core_nt_blast.sh \
  /path/to/barcodes.csv /path/to/new_run_dir /path/to/core_nt/core_nt

.venv/bin/python scripts/audit_complete_blast.py --out /path/to/audited_figure2 \
  --panel B24=/path/to/new_run_dir/raw/B24_blast_raw.tsv \
  --panel K40=/path/to/new_run_dir/raw/K40_blast_raw.tsv \
  --panel K46=/path/to/new_run_dir/raw/K46_blast_raw.tsv \
  --panel K68=/path/to/new_run_dir/raw/K68_blast_raw.tsv

.venv/bin/python scripts/crosscheck_construct_loci.py --audit-dir /path/to/audited_figure2 \
  --barcode-csv /path/to/barcodes.csv --out /path/to/audited_figure2/crosscheck
```

The driver requires BLAST+ 2.15.0 and runs `blastn -task blastn-short -dust no -soft_masking false -max_target_seqs 5000 -evalue 10`. BLAST+ 2.15.0 defaults to an e-value of 1000 for `blastn-short` (2.12.0: 10), so the threshold is passed explicitly. Optional environment variables:

- `PANELS`: a subset of the panels.
- `BARCODES`: comma-separated barcodes for pilot runs.
- `EXPECTED_DB_INFO`: text that must appear in `blastdbcmd -info`, to pin the database release.
- `PYTHON`: the Python interpreter to use.

Each run writes queries, raw TSVs, `db_info.txt` and `run_manifest.json` (versions, commands, hashes) to a new output directory, and refuses to reuse a directory.

Use one complete database release. All volumes should share one date and their summed counts should equal the database metadata (`blastdbcheck` and `blastdbcmd -info` per volume). A partially updated database mixes releases, duplicates records and breaks accession look-ups.

The audit requires `pident = 100`, `mismatch = 0`, `gapopen = 0`, query endpoints spanning 1 to `qlen`, alignment length equal to `qlen`, and an ungapped subject span of the same length. Query length and subject coordinates must be positive. `qcovs` is subject-level coverage and is not proof of individual-HSP completeness. Raw exports use the 17-column schema defined in the audit script.

The audit also checks three things, and exits with status 1 if any fails:

- Forward and reverse-complement queries give identical loci.
- Every construct locus is nested in the shorter panels (K68 ⊆ K46 ⊆ K40 ⊆ B24).
- No query reaches 4,500 subjects, 90 % of the cap.

`construct_context_pairs.csv` lists, for each barcode-only pair, how many of its loci each construct also matches.

The cross-check fetches ±40 nt around every barcode-only locus from NCBI EFetch by versioned accession, rate-limited and cached. It requires the sequence-based construct calls to agree with the BLAST calls. Only accessions and coordinates are sent to NCBI.

Distinct barcode-accession pairs are primary; HSPs are secondary. Grouping retains exported TaxIDs rather than parsing names from record titles. Record-title categories do not independently establish the function or biological origin of a matched interval. The audit writes detailed local HSP/accession working tables as well as aggregate plot tables; keep the detailed outputs out of Git.

`verify_read_constructs.py` checks the constructs in untrimmed reads. It streams the first reads of each FASTQ and reports positional base composition around exact barcode matches as aggregate tables only:

```bash
.venv/bin/python scripts/verify_read_constructs.py --barcode-csv /path/to/barcodes.csv \
  --dataset LABEL=/path/to/untrimmed.fastq --reads 200000 --out /path/to/construct_verification
```

### Taxonomy and GraphPad-Style Figure 2

Map exported TaxIDs to NCBI Taxonomy scientific names, ranks and lineages, then render:

```bash
.venv/bin/python scripts/map_figure2_taxonomy.py \
  --input /path/to/audited_figure2/Fig2A_taxid_groups.csv \
  --out /path/to/audited_figure2 --cache /path/to/taxonomy_cache

Rscript scripts/plot_audited_figure2_graphpad.R \
  --data /path/to/audited_figure2 --out /path/to/rendered_figure2 --tag figure2
```

The first mapping call sends only the TaxIDs to NCBI. Subsequent calls reuse the exact cached XML and request metadata when available, checking its checksum. Missing/duplicate mappings fail explicitly; aliases are recorded without merging the original exported groups. Names and ranks are current NCBI annotations, not independent biological validation or a reconstruction of taxonomy at the BLAST database-build date. Informational species-ancestor columns are not used for regrouping.

The R renderer uses GraphPad-style axes and typography, grey bubbles in panel A and a yellow/teal palette in panel B. It draws panels A and B:

- A: barcode-only TaxID groups, with named and ranked labels and connector lines.
- B: per-barcode accession counts for barcodes with more than three accessions, coloured by the longest native construct matched at the same locus.

Panel C, the title-category composition of the barcode-only pairs of each hit-bearing barcode, is drawn in GraphPad Prism from `Prism_Fig2C_B24_pair_counts.csv` or `Prism_Fig2C_B24_pair_percent.csv`, which `audit_complete_blast.py` writes. The composite leaves its slot empty.

It also writes `Prism_Fig2B_longest_construct.csv` for GraphPad and lists the barcodes with a single accession in `Fig2B_single_accession_barcodes.csv`; panel B omits all barcodes with one to three accessions. It is a study-specific reproduction script: expected totals and selected label positions are deliberately checked. For another dataset, use the general audit/mapping outputs and adapt the figure assertions and layout deliberately.

## Read-Length Histograms

`make_s1_length_histogram.py` exports pooled read-length histograms per condition (by default 10 bp bins from 0 to 2,000 bp, each divided by the largest bin of its condition) from the exact read lengths in each `comparison.sqlite`:

```bash
.venv/bin/python scripts/make_s1_length_histogram.py \
  --cmp-dirs /path/to/experiment_a_demux/cmp_all /path/to/experiment_b_demux/cmp_all \
  --out /path/to/figureS1_histograms
```

Supply every experiment's comparison folder, not just the two illustrative paths. It requires the exact lengths to reproduce each experiment's `length_distribution_bins.csv`. Its `length_counts/` cache holds per-experiment counts and is a local working file, not publication data.

`make_s1_prism_tables.py` defines the coarse length bins used for that reconciliation and exports pooled coarse-bin counts and yield summaries from the per-experiment `length_distribution_bins.csv` files (`--bins`, `--yield-summary`, `--out`). Unequal-width bins need density for a continuous histogram; the open-ended overflow bin has no defined width.

## Adaptive-Sampling Rejections

`summarize_adaptive_sampling.py` counts, per run, the reads that passed MinKNOW's read filter and those among them that ended with an adaptive-sampling rejection (`end_reason` `data_service_unblock_mux_change`), with the median read length of both groups. For all pass reads it also reports the median and arithmetic mean of `mean_qscore_template` and the mean as NanoStat computes it (Q-scores truncated to integers, averaged as error probabilities and converted back). It streams the MinKNOW sequencing summaries line by line:

```bash
.venv/bin/python scripts/summarize_adaptive_sampling.py \
  --summary "experiment_a=/path/to/sequencing_summary_a.txt" \
  --summary "experiment_b=/path/to/sequencing_summary_b.txt" \
  --out /path/to/adaptive_sampling_pass_reads.csv
```

Sequencing summaries contain read identifiers and are local working files, not publication data.

## Residual Barcode Detector

The reported rates use `compare_dorado_barbell_outputs.py --detector occurrences --max-edits 2`. Each barcode orientation is compiled with the Python `regex` package as `(SEQUENCE){e<=2}`, a Levenshtein distance of at most two: any combination of up to two substitutions, insertions and deletions. Without `--max-edits`, the budget is at most one substitution, one insertion and one deletion (`{s<=1,i<=1,d<=1}`).

`OccurrenceScanner` in `validate_detector_occurrences.py` enumerates every candidate occurrence. With `--max-edits 2`, each 24 nt barcode is divided into three 8 nt seeds; two edits can disrupt at most two seeds, so every admissible match contains an exact seed. Seed hits propose starts, shifted by up to 2 nt in either direction, and each 22-26 nt span is confirmed with `fullmatch`. Overlapping candidates across all patterns are resolved greedily by edit distance, deviation from 24 nt, start, end, barcode number and orientation; adjacent matches are kept. `residual_count` is the number of non-overlapping occurrences. Matches within the 80 bp terminal window of either read end are terminal and all others internal; the per-read location classes use the priority internal > both-ends > terminal.

`--detector best_hit` instead keeps one `pattern.search` per barcode orientation (regex BESTMATCH), so it cannot count repeated copies of one orientation; it was not used for the reported rates. Each output directory records `detector_provenance.json` (detector, scanner SHA-256, regex version, edit budget, terminal window), and the scanner refuses to resume into a database made with different settings.

`--decoys N` also scans every read against N decoy barcode sets, in which the bases of each forward barcode are shuffled (composition preserved) and the backward decoy is the reverse complement of the shuffled forward. Decoy matches are located and classified like real ones and summarised in `decoy_background.csv`. This multiplies scan time; `--decoy-seed` fixes the shuffles.

The command line of `validate_detector_occurrences.py` (`--help`) samples reads and synthetic sequences for detector checks; the reported results do not use it.

## Tests

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python -m unittest discover -s scripts/tests -v
```

The 48 tests cover strict HSP criteria, construct-panel construction (sequences, coordinates, nesting, absence of the cDNA primer), the audit's nesting and symmetry checks, taxonomy parsing/alias and failure handling, the occurrence detector (seed candidates, overlap policy and window boundaries; with `--max-edits 2`, candidates equal to an exhaustive dynamic-programming edit-distance comparison, all double deletions and insertions found, three substitutions rejected), synthetic reads that pin the best-hit and occurrence results of the production scanner, the decoy sets, and the adaptive-sampling summary of a synthetic sequencing summary. They do not require study reads, a BLAST database or network taxonomy access. Two tests run only when a 96-barcode `barcodes.csv` is placed in the repository root and are skipped otherwise: the construct check of all 96 barcodes, and the planted-copy check of the decoy scan, which plants every decoy of both sets (seed 20260929) exactly and with one or two edits between random flanks and requires all 3,840 copies to be counted as internal decoy matches with `--max-edits 2`.

## Measurement Limits and Release Notes

- The reported edit budget is a Levenshtein distance of at most two per barcode (`--max-edits 2`). The occurrence detector enumerates non-overlapping occurrences under an operational overlap rule; it does not establish physical barcode copies. One short-read match can satisfy both terminal-window tests. Decoy barcodes estimate coincidental matching within the edit budget but are not matched biological barcode-free controls.
- Condition-level record ratios, barcode-ID multiplicity and location summaries are detector-defined endpoints, not calibrated sample-crosstalk or clinical false-positive rates.
- The legacy `write_stat_tests` branch in the scanner queries `dorado`/`barbell` labels. Its output is not a valid pairwise-test table for the c0-c6 design. Use the separate paired experiment-level analysis; the legacy branch was not silently redefined for the published comparison.
- `stat_mean` length-bin summaries are unweighted experiment means. Pooled percentages require weighting by each condition's output read counts. Geometric summaries are not additive.
- Barcode-like sequences in experimental reads and matching kit barcode sequences in references are complementary observations. They do not establish that the residual reads caused false taxonomic assignments or generated the surveyed reference records.
- NCBI name/rank annotation does not turn mixed-rank bacterial, eukaryotic and viral groups into a bacteria-only species census; title-keyword categories do not establish interval function.

This is a pre-release snapshot. The study's read-level results were produced with Barbell's default pattern set and `--max-edits 2`, and its barcode-database screen with these scripts against one complete core_nt release. The tested environment is not asserted to reproduce the exact historical production environment.

For a reproducible citation, record the Git commit used and cite the Zenodo **version DOI** after the matching release has been archived. A concept DOI does not pin a specific release.