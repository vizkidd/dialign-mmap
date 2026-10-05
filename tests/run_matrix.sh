#!/bin/bash
# run_matrix.sh <orig-binary> <new-binary> [cases-file]
# every line of cases.txt (dataset<TAB>options) is run with both programs and
# all output files are compared byte for byte.
here=$(dirname "$(readlink -f "$0")"); D=$here/data; C=$here/compare.sh; O=$1; N=$2
cf=${3:-$here/cases.txt}; fail=0; n=0
while IFS=$'\t' read -r f opts; do
  [ -z "$f" ] && continue; case $f in \#*) continue;; esac
  n=$((n+1)); "$C" "$O" "$N" "$D/$f.fa" $opts || fail=$((fail+1))
done < "$cf"
echo; [ $fail = 0 ] && echo "ALL $n CASES IDENTICAL" || echo "$fail of $n FAILED"; exit $fail
