#!/bin/bash
# robustness.sh <new-binary-source-dir> [cases-file]
# Detects undefined behaviour / uninitialised reads / use-after-free that the
# plain comparison against the original can miss.  The new program is built in
# several ways and every case must give the SAME output in all of them:
#   default | -nommap + MALLOC_PERTURB_ | DIALIGN_MMDEBUG (poisoned frees) |
#   -ftrivial-auto-var-init=zero | =pattern | -O0 | -threads 3
here=$(dirname "$(readlink -f "$0")"); src=$(readlink -f "$1"); cf=${2:-$here/cases.txt}
export DIALIGN2_DIR=$here/../dialign2_dir
B=$(mktemp -d)
b(){ gcc -O2 -g -std=gnu99 -ffp-contract=off -DCONS -fopenmp $2 -I"$src" "$src"/*.c -o "$B/$1" -lm -lpthread 2>/dev/null; }
b base ""; b zero "-ftrivial-auto-var-init=zero"; b pattern "-ftrivial-auto-var-init=pattern"; b o0 "-O0"
fail=0; n=0
while IFS=$'\t' read -r f opts; do
  [ -z "$f" ] && continue; case $f in \#*) continue;; esac
  case "$opts" in *-ts*|*-stdo*) continue;; esac     # -ts prints wall-clock times
  n=$((n+1)); sums=""
  for v in "base|" "base|-nommap|165" "base||dbg" "zero|" "pattern|" "o0|" "base|-threads 3"; do
    IFS='|' read -r bin extra pert <<< "$v"
    w=$(mktemp -d); cp "$here/data/$f.fa" "$w/"; for x in anc xfr; do cp "$here/data/$f.$x" "$w/" 2>/dev/null; done
    cp "$B/$bin" "$w/d"
    ( cd "$w"; if [ "$pert" = dbg ]; then export DIALIGN_MMDEBUG=1; elif [ -n "$pert" ]; then export MALLOC_PERTURB_=$pert; fi
      ./d $opts $extra $f.fa >stdout 2>/dev/null; echo rc=$? >> stdout )
    sums="$sums $( (cd $w; for x in $(ls | grep -v '^d$\|\.fa$\|\.anc$\|\.xfr$'); do echo "== $x"; grep -v 'program call\|program parameters\|dialign2\|\./d ' $x; done) | md5sum | cut -c1-8)"
    rm -rf "$w"
  done
  u=$(echo $sums | tr ' ' '\n' | sort -u | wc -l)
  if [ "$u" != 1 ]; then echo "UNSTABLE: $f [$opts] :$sums"; fail=$((fail+1)); fi
done < "$cf"
rm -rf "$B"; echo; [ $fail = 0 ] && echo "ROBUST: $n cases give identical output in all 7 build/run variants" || echo "$fail of $n UNSTABLE"; exit $fail
