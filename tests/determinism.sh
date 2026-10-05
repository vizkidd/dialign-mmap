#!/bin/bash
# determinism.sh <orig-binary> <new-binary>
# the new program must give the same result for 1..8 threads, with and
# without the mmap store, with a different temp dir, and on repeated runs
here=$(dirname "$(readlink -f "$0")"); D=$here/data; C=$here/compare.sh; O=$(readlink -f "$1"); N=$(readlink -f "$2"); fail=0
tmp=$(mktemp -d)
for th in 1 2 3 4 8; do
  for f in prot40 dna36 trna_50 ties30 trnaS_100 trnaL_40; do
    OMP_NUM_THREADS=$th "$C" "$O" "$N" "$D/$f.fa" -n >/dev/null || { echo "FAIL threads=$th $f"; fail=$((fail+1)); }
  done
done
for f in trna_50 ties30 dna36; do
  for th in 2 3 8; do OMP_NUM_THREADS=$th "$C" "$O" "$N" "$D/$f.fa" -n -ow >/dev/null || { echo "FAIL -ow threads=$th $f"; fail=$((fail+1)); }; done
done
for opt in "-nommap" "-tmpdir $tmp" "-threads 3" "-nommap -threads 2"; do
  # (the original rejects the new options, so compare new-vs-new instead)
  export DIALIGN2_DIR=$here/../dialign2_dir
  a=$(mktemp -d); b=$(mktemp -d); cp "$D/prot15.fa" $a; cp "$D/prot15.fa" $b
  (cd $a && "$N" prot15.fa >/dev/null); (cd $b && "$N" $opt prot15.fa >/dev/null)
  [ -s $a/prot15.ali ] && [ -s $b/prot15.ali ] || { echo "FAIL [$opt]: no output"; fail=$((fail+1)); rm -rf $a $b; continue; }
  diff <(grep -v "program call" $a/prot15.ali) <(grep -v "program call" $b/prot15.ali) >/dev/null || { echo "FAIL [$opt]"; fail=$((fail+1)); }
  rm -rf $a $b
done
[ -z "$(ls -A $tmp)" ] || { echo "FAIL: temp dir not clean"; fail=$((fail+1)); }
rm -rf $tmp
[ $fail = 0 ] && echo "DETERMINISM OK" || echo "$fail FAILED"; exit $fail
