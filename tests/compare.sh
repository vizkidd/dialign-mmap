#!/bin/bash
# compare.sh <orig-binary> <new-binary> <file.fa> [dialign options...]
#
# Runs both programs on a private copy of the input (same argv[0] and same
# arguments, so even the "program call:" line matches) and compares every
# file the ORIGINAL program wrote - the .ali alignment, its stdout and all
# optional outputs (.fa .msf .cw .frg .fop .afc .fsm .log .mat .cap ...).
# Exit status 0 = byte-identical.
ORIG=$(readlink -f "$1"); NEW=$(readlink -f "$2"); f=$(readlink -f "$3"); shift 3
export DIALIGN2_DIR=${DIALIGN2_DIR:-$(dirname "$(readlink -f "$0")")/../dialign2_dir}
b=$(basename "$f" .fa); src=$(dirname "$f")
w=$(mktemp -d); mkdir "$w/o" "$w/n"
for d in o n; do
  cp "$f" "$w/$d/"; for x in anc xfr; do [ -f "$src/$b.$x" ] && cp "$src/$b.$x" "$w/$d/"; done
done
cp "$ORIG" "$w/o/dialign2-2"; cp "$NEW" "$w/n/dialign2-2"
s=$(date +%s.%N); (cd "$w/o" && ./dialign2-2 "$@" "$b.fa" >stdout 2>stderr); ro=$?
m=$(date +%s.%N); (cd "$w/n" && ./dialign2-2 "$@" "$b.fa" >stdout 2>stderr); rn=$?
e=$(date +%s.%N)
status=IDENTICAL
for x in "$w"/o/*; do
  n=$(basename "$x"); case $n in stderr|dialign2-2) continue;; esac
  cmp -s "$x" "$w/n/$n" || { status="DIFF($n)"; }
done
# the new program must not leave temporary files behind
for x in "$w"/n/*; do n=$(basename "$x"); [ -e "$w/o/$n" ] || [ "$n" = stderr ] || status="$status EXTRA($n)"; done
if [ "$ro" -gt 128 ] && [ "$rn" = 0 ]; then status="ORIG-CRASHES(rc=$ro);NEW-COMPLETES"; fi
# expected, intentional divergence: KNOWN_ISSUES.md #1 fix. The original
# silently clobbers a "*.fa"-named input file when run with -fa; the port
# refuses to and writes the alignment to "<input>.dialign-aligned.fa"
# instead, so the only "difference" is that the input survives.
if [ "$status" = "DIFF($b.fa) EXTRA($b.fa.dialign-aligned.fa)" ]; then
  status="ORIG-CLOBBERS-INPUT(fixed, see KNOWN_ISSUES.md);NEW-PRESERVES-INPUT"
fi
printf "%-16s %-34s orig rc=%s %6.1fs | new rc=%s %6.1fs | %s\n" "$b" "$*" $ro $(echo "$m-$s"|bc) $rn $(echo "$e-$m"|bc) "$status"
rm -rf "$w"
case "$status" in ORIG-CRASHES*|ORIG-CLOBBERS-INPUT*) exit 0;; esac
[ "$status" = IDENTICAL ] && [ $ro = $rn ]
