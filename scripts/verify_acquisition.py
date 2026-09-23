"""Compare imported GTEx all-pairs values with the original public API records."""
import csv
import gzip
import json
import math
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
with open(ROOT / 'results/public/variant_gene_evidence.tsv', encoding='utf-8', newline='') as handle:
    evidence = list(csv.DictReader(handle, delimiter='\t'))
expected = {(r['gene_id'].split('.')[0], r['tissue'], r['variant_id'].removesuffix('_b38')): r for r in evidence}
hits = set()
for path in (ROOT / 'data/raw/summary_statistics/regions').glob('Colon_*/*.tsv.gz'):
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = (row['gene_id'].split('.')[0], path.parent.name, row['variant'])
            if key not in expected:
                continue
            ref = expected[key]
            if not math.isclose(float(row['pvalue']), float(ref['eqtl_p']), rel_tol=1e-4):
                raise ValueError(f'eQTL p-value differs from original API: {key}')
            if not math.isclose(float(row['beta']), float(ref['nes']), rel_tol=1e-4, abs_tol=1e-8):
                raise ValueError(f'eQTL effect differs from original API: {key}')
            hits.add(key)
missing = sorted(set(expected) - hits)
report = dict(verified_at_utc=datetime.now(timezone.utc).isoformat(),
              expected_associations=len(expected), matched_associations=len(hits),
              missing=missing, relative_tolerance=1e-4, beta_absolute_tolerance=1e-8,
              status='passed' if not missing else 'incomplete')
output = ROOT / 'data/processed/acquisition/original_api_validation.json'
output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(json.dumps(report, indent=2))
if missing:
    raise SystemExit('Acquired data do not contain every original API association')
