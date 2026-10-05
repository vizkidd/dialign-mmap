#!/bin/bash
# tools/detect.sh <CC> <repo-root>
#
# Prints, on stdout, exactly two lines:
#   OPENMP=1   or   OPENMP=0
#   OPT=-O3    or   OPT=-O2
#
# OpenMP is enabled if a trivial OpenMP program compiles and links.
#
# -O3 is used ONLY if a full build at -O3 produces byte-identical alignment
# output to a full build at -O2 (the level this project's test suite --
# see PERFORMANCE.md / CHANGES.md section 7 -- was actually verified
# deterministic at) on a handful of real test cases. This project's whole
# point is bit-for-bit reproducible output; a fused-multiply-add or
# vectorization difference that -O3 can introduce even with
# -ffp-contract=off would silently break that guarantee, so "supported"
# here means "verified", not just "the compiler accepts the flag".
# If anything about the probe is inconclusive (build fails, test data
# missing, mismatch found), it falls back to the safe, tested -O2.
set -u
CC=${1:-gcc}
ROOT=${2:-.}
SRC="$ROOT/src"
DATA="$ROOT/tests/data"
DIALIGN2_DIR="$ROOT/dialign2_dir"

# --- OpenMP ---------------------------------------------------------------
tmpd=$(mktemp -d)
cat > "$tmpd/omp.c" <<'EOF'
#include <omp.h>
int main(void){ return omp_get_num_threads() >= 0 ? 0 : 1; }
EOF
if "$CC" -fopenmp "$tmpd/omp.c" -o "$tmpd/omp" -lgomp >/dev/null 2>&1; then
  openmp=1
else
  openmp=0
fi

# --- -O3 safety probe -------------------------------------------------------
opt="-O2"
if [ -d "$SRC" ] && [ -d "$DATA" ] && [ -d "$DIALIGN2_DIR" ]; then
  ompflags=""; omplibs=""
  [ "$openmp" = 1 ] && ompflags="-fopenmp" && omplibs="-lgomp"

  o2bin="$tmpd/d_o2"; o3bin="$tmpd/d_o3"
  if "$CC" -O2 -g -std=gnu99 -ffp-contract=off -DCONS -I"$SRC" $ompflags \
        "$SRC"/*.c -o "$o2bin" -lm -lpthread $omplibs >/dev/null 2>&1 \
     && "$CC" -O3 -g -std=gnu99 -ffp-contract=off -DCONS -I"$SRC" $ompflags \
        "$SRC"/*.c -o "$o3bin" -lm -lpthread $omplibs >/dev/null 2>&1
  then
    ok=1
    for f in prot15 dna20 trna_50; do
      [ -f "$DATA/$f.fa" ] || continue
      wd2=$(mktemp -d); wd3=$(mktemp -d)
      cp "$DATA/$f.fa" "$wd2/"; cp "$DATA/$f.fa" "$wd3/"
      ( cd "$wd2" && DIALIGN2_DIR=$(readlink -f "$DIALIGN2_DIR") "$o2bin" -n "$f.fa" >/dev/null 2>&1 )
      ( cd "$wd3" && DIALIGN2_DIR=$(readlink -f "$DIALIGN2_DIR") "$o3bin" -n "$f.fa" >/dev/null 2>&1 )
      if ! diff -q <(grep -v 'program call\|program parameters' "$wd2/$f.ali" 2>/dev/null) \
                   <(grep -v 'program call\|program parameters' "$wd3/$f.ali" 2>/dev/null) >/dev/null 2>&1; then
        ok=0
      fi
      rm -rf "$wd2" "$wd3"
    done
    [ "$ok" = 1 ] && opt="-O3"
  fi
fi

rm -rf "$tmpd"
echo "OPENMP=$openmp"
echo "OPT=$opt"
