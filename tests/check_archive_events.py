"""Independent server-side request count and timing oracles."""
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
record = json.loads((root/'events.json').read_text())
counts = record['counts']
assert counts['retry429'] == 2, counts
assert counts['date429'] == 2, counts
assert counts['retry500'] == 3, counts
assert counts['drop'] == 2, counts
assert counts['missing'] == 1, counts
assert counts['oversize'] == 2, counts
assert counts['truncated'] == 1, counts
assert counts['chunked'] == 1, counts
assert counts['slow'] == 1, counts
assert counts['split'] == 6, counts  # two uncached reads; the cache hit sent none
assert counts['always500'] == 4, counts
assert counts['windowretry'] == 3, counts
versioned = [e['url'] for e in record['events'] if e['mode'] == 'versioned']
assert len(versioned) == 3, versioned
assert 'version=latest' in versioned[0], versioned
assert all('version=v1' in url for url in versioned[1:]), versioned
backoff = [e['time'] for e in record['events'] if e['mode'] == 'retry500']
assert backoff[1] - backoff[0] >= 0.03, backoff
assert backoff[2] - backoff[1] >= 0.06, backoff
assert 'deadline' not in counts, counts
parallel = [e['time'] for e in record['events'] if e['mode'] == 'parallel']
assert len(parallel) == 4, parallel
separation = min(b-a for a,b in zip(parallel,parallel[1:]))
assert separation >= 0.1, separation
print(f'parallel_min_start_separation_s={separation:.6f}')
assert all('TIAGO/' in e['agent'] for e in record['events'])
for mode, floor in [('retry429', 1.0), ('date429', 1.0), ('interval', 0.5)]:
    t = [e['time'] for e in record['events'] if e['mode'] == mode]
    delta = t[1]-t[0]
    assert delta >= floor, (mode, delta)
    print(f'{mode}_server_interval_s={delta:.6f}')
print('Independent request-count, User-Agent and timing checks passed')
