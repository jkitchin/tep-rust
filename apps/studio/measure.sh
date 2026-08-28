#!/usr/bin/env bash
#
# Report what TEP Studio costs to load, raw and gzipped.
#
# PLAN.org budgets the whole app at under 1.5 MB gzipped. A budget nobody
# measures is a wish, and the thing that eats it is never the wasm: it is one
# convenient JavaScript dependency. So this prints every file that the browser
# has to fetch before the first frame, gzipped, because gzip is what a static
# host actually serves.
#
# Usage: apps/studio/measure.sh [dist-directory]

set -euo pipefail

dist="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/dist}"

if [ ! -d "${dist}" ]; then
  echo "no such directory: ${dist}" >&2
  echo "run apps/studio/build.sh first" >&2
  exit 1
fi

total_raw=0
total_gz=0

printf '%-34s %12s %12s\n' 'file' 'raw' 'gzip -9'
printf '%-34s %12s %12s\n' '----' '---' '-------'

# Sorted so the report is stable across machines and comparable between runs.
while IFS= read -r file; do
  rel="${file#"${dist}"/}"
  raw=$(wc -c < "${file}" | tr -d ' ')
  gz=$(gzip -9 -c "${file}" | wc -c | tr -d ' ')
  total_raw=$((total_raw + raw))
  total_gz=$((total_gz + gz))
  printf '%-34s %12s %12s\n' "${rel}" "${raw}" "${gz}"
done < <(find "${dist}" -type f \
  \( -name '*.wasm' -o -name '*.js' -o -name '*.html' -o -name '*.css' \) | sort)

printf '%-34s %12s %12s\n' '----' '---' '-------'
printf '%-34s %12s %12s\n' 'total' "${total_raw}" "${total_gz}"

budget=$((1500 * 1024))
printf '\nbudget %s bytes gzipped (PLAN.org: under 1.5 MB): ' "${budget}"
if [ "${total_gz}" -le "${budget}" ]; then
  printf 'met, %s%% used\n' "$((total_gz * 100 / budget))"
else
  printf 'EXCEEDED by %s bytes\n' "$((total_gz - budget))"
  exit 1
fi
