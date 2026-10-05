#!/bin/bash
# tree_unit.sh <original-functions.c>
# Compares the original av_tree_print() with the rewritten one on random,
# tie-heavy similarity matrices (all three linkage modes, N up to 250).
here=$(dirname "$(readlink -f "$0")"); src=$here/../src; w=$(mktemp -d)
python3 - "$1" "$src/functions.c" "$w" <<'PY'
import sys
def ext(path,name):
    s=open(path).read(); a=s.index('void av_tree_print'); b=s.index('void print_log(',a)
    return s[a:b].replace('av_tree_print','av_tree_print_'+name,1)
open(sys.argv[3]+'/orig_fn.inc','w').write(ext(sys.argv[1],'orig'))
open(sys.argv[3]+'/new_fn.inc','w').write(ext(sys.argv[2],'new'))
PY
cp "$here/tree_unit.c" "$w/t.c"
gcc -O1 -std=gnu99 -I"$src" -I"$w" "$w/t.c" "$src/mmstore.c" -o "$w/t" -lm -lpthread || exit 2
rc=0; for s in 1 2 3; do "$w/t" $s || rc=1; done; rm -rf "$w"; exit $rc
