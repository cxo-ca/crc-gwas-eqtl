"""Durable Windows download job; reuses the validated acquisition pipeline."""
import json
import msvcrt
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT / 'data/processed/acquisition/download_job_status.json'


def now():
    return datetime.now(timezone.utc).isoformat()


def update(state):
    state['updated_at_utc'] = now()
    state['completed_gwas_windows'] = len(list((ROOT / 'data/raw/summary_statistics/regions/GCST90129505').glob('*.tsv.gz')))
    state['completed_eqtl_windows'] = len(list((ROOT / 'data/raw/summary_statistics/regions').glob('Colon_*/*.tsv.gz')))
    state['planned_eqtl_windows'] = 61
    temporary = STATUS.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(state, indent=2) + '\n', encoding='utf-8')
    temporary.replace(STATUS)


def main():
    os.chdir(ROOT)
    STATUS.parent.mkdir(parents=True, exist_ok=True)
    lock = open(STATUS.with_suffix('.lock'), 'a+b')
    lock.seek(0)
    if not lock.read(1):
        lock.write(b'0')
        lock.flush()
    lock.seek(0)
    msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
    state = dict(status='running', pid=os.getpid(), started_at_utc=now(),
                 max_attempts=12, timeout_hours=8)
    deadline = time.monotonic() + 8 * 3600
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    for attempt in range(1, 13):
        log = ROOT / 'logs' / f'download_job_{stamp}_{attempt:02d}.log'
        state.update(attempt=attempt, log=log.relative_to(ROOT).as_posix(), status='running')
        update(state)
        with open(log, 'wb') as output:
            child = subprocess.Popen([sys.executable, '-u', 'scripts/acquire_regional_sumstats.py'],
                                     cwd=ROOT, stdout=output, stderr=subprocess.STDOUT)
            state['worker_pid'] = child.pid
            update(state)
            while child.poll() is None:
                if time.monotonic() >= deadline:
                    child.terminate()
                    child.wait()
                    state.update(status='failed', reason='Eight-hour job time limit reached; cache retained')
                    update(state)
                    return 1
                time.sleep(10)
                update(state)
        state['exit_code'] = child.returncode
        if child.returncode == 0:
            summary_path = STATUS.parent / 'acquisition_summary.json'
            summary = json.loads(summary_path.read_text(encoding='utf-8'))
            state.update(status='completed_available_data', summary=summary_path.relative_to(ROOT).as_posix(),
                         paired_genes=summary['paired_genes'], unavailable_gwas_regions=summary['unavailable_gwas_regions'])
            update(state)
            return 0
        message = log.read_text(encoding='utf-8', errors='replace')
        if not re.search(r'curl range request failed \((6|7|28|35|52|56)\)', message):
            state.update(status='failed', reason='Non-network failure: inspect log; not automatically ignored')
            update(state)
            return child.returncode
        state.update(status='retrying_network', reason='Transient transfer failure; validated cache retained')
        update(state)
        time.sleep(30)
    state.update(status='failed', reason='Network retry limit reached; validated cache retained')
    update(state)
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
