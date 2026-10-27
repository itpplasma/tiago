"""Loopback-only CTest fixture. Every signal and configuration is synthetic."""
import datetime
import email.utils
import http.server
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time
import urllib.parse
import urllib.request


def serve(root):
    lock = threading.Lock()
    counts = {}
    window_counts = {}
    events = []
    active = 0
    peak = 0

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *_):
            pass

        def do_GET(self):
            nonlocal active, peak
            url = urllib.parse.urlsplit(self.path)
            mode = url.path.lstrip('/')
            if mode == 'stats':
                with lock:
                    body = json.dumps(dict(counts=counts, events=events, peak=peak)).encode()
                self.send_response(200)
                self.end_headers()
                self.wfile.write(body)
                return
            with lock:
                counts[mode] = counts.get(mode, 0) + 1
                count = counts[mode]
                active += 1
                peak = max(peak, active)
                events.append(dict(mode=mode, time=time.monotonic(),
                                   agent=self.headers.get('User-Agent', ''), url=self.path))
                (root / 'events.json').write_text(json.dumps(dict(counts=counts, events=events, peak=peak)))
            try:
                args = urllib.parse.parse_qs(url.query)
                first = float(args.get('first', ['0'])[0])
                last = float(args.get('last', ['1'])[0])
                version = args.get('version', ['v1'])[0]
                if mode == 'windowretry':
                    key = (first, last)
                    with lock:
                        window_counts[key] = window_counts.get(key, 0) + 1
                        attempt = window_counts[key]
                    if attempt == 1:
                        self.send_response(503)
                        self.end_headers()
                        return
                if mode == 'versioned' and version == 'latest':
                    version = 'v1' if count == 1 else 'v2'
                if mode in ('retry429', 'date429') and count == 1:
                    self.send_response(429)
                    retry = '1'
                    if mode == 'date429':
                        retry = email.utils.format_datetime(datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(seconds=2), usegmt=True)
                    self.send_header('Retry-After', retry)
                    self.end_headers()
                    return
                if mode == 'always500' or (mode == 'retry500' and count < 3):
                    self.send_response(503)
                    self.end_headers()
                    return
                if mode == 'drop' and count == 1:
                    self.close_connection = True
                    return
                if mode == 'missing':
                    self.send_response(404)
                    self.end_headers()
                    return
                if mode in ('slow', 'parallel'):
                    time.sleep(0.3 if mode == 'slow' else 0.1)
                n = 20000 if mode == 'large' else 3
                times = [first + (last-first)*i/(n-1) for i in range(n)]
                values = [2*t+1 for t in times]
                result = dict(data=dict(time=times, values=values, variance=[0.25]*n), version=version)
                if mode == 'wrongversion':
                    result['version'] = 'unexpected'
                if mode == 'malformed':
                    result['data']['values'] = [1.0]
                if mode == 'outside':
                    result['data']['time'][0] = first - 1
                if mode == 'null':
                    result['data']['values'][0] = None
                if mode == 'negative':
                    result['data']['variance'][0] = -1.0
                if mode == 'empty':
                    result['data'] = dict(time=[], values=[], variance=[])
                body = json.dumps(result).encode()
                if mode == 'oversize':
                    body = b'x'*65536
                if mode == 'badjson':
                    body = b'{"broken":'
                if mode == 'binary':
                    body = bytes([0, 1, 127, 128, 255, 0])
                self.send_response(200)
                self.send_header('Content-Type', 'application/octet-stream' if mode == 'binary' else 'application/json')
                if mode == 'chunked':
                    self.send_header('Transfer-Encoding', 'chunked')
                if mode == 'truncated':
                    self.send_header('Content-Length', str(len(body)))
                    body = body[:len(body)//2]
                # Oversize has no length: enforce the cap from streamed bytes.
                self.end_headers()
                for i in range(0, len(body), 73):
                    chunk = body[i:i+73]
                    if mode == 'chunked':
                        self.wfile.write(f'{len(chunk):x}\r\n'.encode()+chunk+b'\r\n')
                    else:
                        self.wfile.write(chunk)
                if mode == 'chunked':
                    self.wfile.write(b'0\r\n\r\n')
            except (BrokenPipeError, ConnectionResetError):
                pass
            finally:
                with lock:
                    active -= 1

        def do_POST(self):
            if self.path != '/shutdown':
                self.send_error(404)
                return
            self.send_response(200)
            self.end_headers()
            threading.Thread(target=self.server.shutdown).start()

    server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    base = f'http://127.0.0.1:{server.server_port}'
    (root / 'base.txt').write_text(base)
    modes = ('ok', 'split', 'retry429', 'date429', 'retry500', 'always500',
             'oversize', 'slow', 'wrongversion', 'malformed', 'outside', 'null',
             'empty', 'badjson', 'binary', 'drop', 'missing', 'large',
             'windowretry', 'versioned', 'chunked', 'truncated', 'negative')
    channels = [dict(logical=mode, url=base+'/'+mode+'?first={first}&last={last}&version={version}',
                     units='eV', time_path='data.time', value_path='data.values',
                     variance_path='data.variance', version_path='version',
                     version='latest' if mode == 'versioned' else 'v1', scale=2.0, offset=1.0) for mode in modes]
    config = dict(source_id='synthetic-loopback', channels=channels,
                  policy=dict(min_interval=0.01, backoff=0.02, max_backoff=0.2,
                              jitter=0.2, max_window=1.0, timeout=2.0, max_elapsed=10.0))
    (root / 'adapter.json').write_text(json.dumps(config))
    (root / 'ready').write_text('ready')
    server.serve_forever()
    server.server_close()


def main():
    action, directory = sys.argv[1:]
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    if action == 'serve':
        serve(root)
    elif action == 'start':
        (root / 'ready').unlink(missing_ok=True)
        log = (root / 'server.log').open('w')
        proc = subprocess.Popen([sys.executable, __file__, 'serve', str(root)],
                                stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                                start_new_session=True)
        for _ in range(200):
            if (root / 'ready').exists():
                return
            if proc.poll() is not None:
                raise RuntimeError('mock server exited; see server.log')
            time.sleep(0.05)
        proc.terminate()
        raise RuntimeError('mock server did not start')
    elif action == 'stop':
        base = (root / 'base.txt').read_text()
        with urllib.request.urlopen(urllib.request.Request(base+'/shutdown', method='POST'), timeout=5):
            pass
    else:
        raise ValueError(action)


if __name__ == '__main__':
    main()
