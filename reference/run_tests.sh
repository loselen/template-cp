#!/usr/bin/env bash
# Run every reference file's self-tests. Exit status is non-zero if any fail.
cd "$(dirname "$0")" || exit 1
PY=${PYTHON:-python3}
fail=0
for f in [0-9][0-9]_*.py; do
  start=$(date +%s%N)
  if out=$(timeout 300 "$PY" "$f" 2>&1) && grep -q "ALL OK" <<<"$out"; then
    ms=$(( ($(date +%s%N) - start) / 1000000 ))
    printf "PASS %-24s %d.%01ds\n" "$f" $((ms / 1000)) $((ms % 1000 / 100))
  else
    printf "FAIL %s\n%s\n" "$f" "$(tail -n 20 <<<"$out")"
    fail=1
  fi
done
exit $fail
