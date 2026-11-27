#!/usr/bin/env python3
"""Export Supplementary Data 1 from an existing, verified Paper 3 run; no searches are run."""
import argparse
import csv
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path


def read_rows(path, delimiter='\t'):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export(figure_sources, search_status, benchmark_by_length, output, genbank_record=None):
    figure_sources, search_status, benchmark_by_length, output = map(Path, (figure_sources, search_status, benchmark_by_length, output))
    manifest_path = figure_sources / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest['contract'] != 'publication_figure_sources_v1':
        raise ValueError('Expected verified publication figure sources')
    for name in ['logical_queries.tsv', 'position_sensitivity.tsv']:
        if sha(figure_sources / name) != manifest['files'][name]:
            raise ValueError('Figure-source checksum differs: ' + name)
    lengths = {}
    for row in read_rows(figure_sources / 'logical_queries.tsv'):
        identifier, length = row['sequence_id'], int(row['length'])
        if identifier in lengths and lengths[identifier] != length:
            raise ValueError('Inconsistent sequence length: ' + identifier)
        lengths[identifier] = length
    searches = read_rows(search_status)
    if len({r['query_id'] for r in searches}) != len(searches) or {r['query_id'] for r in searches} != set(lengths):
        raise ValueError('Search rows must cover every distinct sequence exactly once')
    for row in searches:
        if row['status'] not in ['success', 'success_zero', 'success_hits']:
            raise ValueError('Incomplete search: ' + row['query_id'])
        row['length_nt'] = length = lengths[row['query_id']]
        row['length_class'] = ('8, 10, or 16 nt' if length in [8, 10, 16] else
                               '19–23 nt' if 19 <= length <= 23 else
                               '24–31 nt' if 24 <= length <= 31 else '≥32 nt')
    positions = read_rows(figure_sources / 'position_sensitivity.tsv')
    for row in positions:
        row['unit_rule'] = row['unit_rule'].replace('primary strict event', 'strict event')
    if output.exists():
        raise FileExistsError('Use a new output directory: ' + str(output))
    output.mkdir(parents=True)
    def write(name, rows):
        with (output / name).open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    write('distinct_searches.csv', searches)
    counts = Counter((r['length_class'], int(r['hsp_rows']) > 0) for r in searches)
    write('distinct_search_counts.csv', [{'length_class': c, 'result': 'Raw HSPs' if hit else 'No HSP',
                                         'searches': counts[(c, hit)]}
          for c in ['8, 10, or 16 nt', '19–23 nt', '24–31 nt', '≥32 nt'] for hit in [True, False]])
    write('positional_sensitivity_all_settings.csv', positions)
    shutil.copyfile(benchmark_by_length, output / 'benchmark_pair_recall_by_length_all_sizes.csv')
    inputs = [manifest_path, figure_sources / 'logical_queries.tsv', figure_sources / 'position_sensitivity.tsv', search_status, benchmark_by_length]
    if genbank_record:
        genbank_record = Path(genbank_record)
        if 'VERSION     AY738638.1' not in genbank_record.read_text() or 'derived from bacteriophage Mu' not in genbank_record.read_text():
            raise ValueError('Expected the checked AY738638.1 GenBank record')
        shutil.copyfile(genbank_record, output / 'AY738638.1.gb'); inputs.append(genbank_record)
    benchmark = read_rows(benchmark_by_length, ',')
    summary = {'distinct_searches': len(searches), 'raw_HSP_searches': sum(int(r['hsp_rows']) > 0 for r in searches),
               'no_HSP_searches': sum(int(r['hsp_rows']) == 0 for r in searches), 'positional_combinations': len(positions),
               'benchmark_database_sizes': sorted({r['database'] for r in benchmark})}
    readme = '''SUPPLEMENTARY DATA 1 — export of the existing release 1.0.1 run

No searches or scientific analyses are rerun by this exporter.

distinct_searches.csv: one row per distinct searched sequence. query_id identifies the sequence, batch_id its BLAST batch, status its completion state, hsp_rows the raw HSP count, unique_targets the distinct subjects, strict_targets/tolerant_targets qualifying subjects, duplicate_hsp_rows removed duplicates, target_limit_reached the search-status flag, length_nt the sequence length, and length_class its disjoint Table S5D class. Multiple panel/control rows may share one search. Raw HSPs need not qualify as matches.

distinct_search_counts.csv: Table S5D totals by length class and raw-HSP result, with each distinct sequence counted once.

positional_sensitivity_all_settings.csv: every combination of platform/global group, window, minimum record length and event-selection variant. Window and minimum_length are in nt. Events and terminal are observed counts, expected is the sum of uniform-placement terminal probabilities, and ratio is observed/expected. A zero-event combination has no estimable ratio. Hash selection uses event identifiers and coordinates without ranking by end distance. units_without_project refers to all eligible record lengths, not each thresholded subset.

benchmark_pair_recall_by_length_all_sizes.csv: recovery by length, arm, rule and database size (core_nt, one tenth and ten times). truth_positive_pairs is the sequence-defined denominator; recovered_pairs and missed_pairs partition it; pair_recall is the recovered fraction. Diagnostic arms use E = 10^12 and calibrated arms E = 10. Table S9B in the Word supplement gives core_nt only.

AY738638.1.gb, when supplied: exact-version GenBank record checked for the pJK95 vector annotation. repeat_region 20..40 is labelled "derived from bacteriophage Mu", within a mini-Mu mobile_element at 20..3310. This supports vector-context interpretation of the overlapping 20..69 match.

provenance.json: input and output SHA-256 checksums plus export totals. See the accompanying manuscript for matching rules, testability and limitations.
'''
    (output / 'README.txt').write_text(readme, encoding='utf-8')
    provenance = {'release': '1.0.1', 'exporter_sha256': sha(__file__), 'summary': summary,
                  'inputs': [{'filename': p.name, 'sha256': sha(p)} for p in inputs],
                  'files': {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()}}
    (output / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--figure-sources', type=Path, required=True)
    p.add_argument('--search-status', type=Path, required=True)
    p.add_argument('--benchmark-by-length', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--genbank-record', type=Path)
    a = p.parse_args()
    print(json.dumps(export(a.figure_sources, a.search_status, a.benchmark_by_length, a.output, a.genbank_record), indent=2))
