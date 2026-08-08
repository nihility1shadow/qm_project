#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source_file="$script_dir/convergence_analyzer.cpp"
binary="$script_dir/.convergence_analyzer.bin"

if ! command -v g++ >/dev/null 2>&1; then
  echo "error: g++ was not found; load the GCC module first" >&2
  exit 2
fi

if [[ ! -x "$binary" || "$source_file" -nt "$binary" ]]; then
  g++ -O3 -std=c++17 -DNDEBUG "$source_file" -o "$binary"
fi

exec "$binary" "$@"
