"""Invalidate mismatching cache entries; original bytes remain for inspection."""
from pathlib import Path
import hashlib
import json

base = Path('data/raw/summary_statistics/range_cache')
report = []
for path in base.glob('*/*/manifest.json'):
    manifest = json.loads(path.read_text(encoding='utf-8'))
    good = {}
    for entry in manifest['ranges']:
        valid = True
        for name_key, hash_key in [('file', 'sha256'), ('headers', 'headers_sha256')]:
            item = path.parent / entry[name_key]
            if not item.exists() or (entry.get(hash_key) and hashlib.sha256(item.read_bytes()).hexdigest() != entry[hash_key]):
                valid = False
        if valid:
            good[entry['file']] = entry
        else:
            report.append({'manifest': str(path), 'file': entry['file'], 'action': 'invalidate_and_redownload'})
    if len(good) != len(manifest['ranges']):
        path.with_suffix('.before_repair.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        manifest['ranges'] = list(good.values())
        path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
Path('logs/range_cache_repair.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print('Invalidated entries:', len(report))
