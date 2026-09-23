import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

folder = Path('data/processed/acquisition')
paths = [folder / 'gwas.regional.tsv.gz', folder / 'eqtl.regional.tsv.gz', folder / 'gtex_expression_sdY.tsv']


def fingerprint(path):
    digest = hashlib.sha256()
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


before = {str(path): fingerprint(path) for path in paths}
with open('logs/offline_acquisition_replay.log', 'wb') as log:
    for script in ['scripts/fetch_gtex_expression.py', 'scripts/acquire_regional_sumstats.py']:
        subprocess.run([sys.executable, script, '--offline'], stdout=log, stderr=subprocess.STDOUT, check=True)
after = {str(path): fingerprint(path) for path in paths}
if before != after:
    raise SystemExit('Offline replay changed the uncompressed scientific data')
report = dict(status='passed', verified_at_utc=datetime.now(timezone.utc).isoformat(),
              comparison='SHA-256 of uncompressed TSV content; gzip timestamps excluded', hashes=after)
(folder / 'offline_validation.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(json.dumps(report, indent=2))
