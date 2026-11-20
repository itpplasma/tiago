#!/bin/sh
# Provider stand-in for fast tests: answers every request with one G-EQDSK.
# Usage: geqdsk_stub_provider.sh GEQDSK REQUEST OUTPUT
set -eu
cp "$1" "$3"
printf 'converged 1\n' > "$3.status"
