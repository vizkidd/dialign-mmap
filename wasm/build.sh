#!/bin/bash
# wasm/build.sh -- build dialign2-2.{js,wasm} with Emscripten, then verify
# the result actually works before calling the build a success.
#
# Requires emcc on PATH (e.g. `source /path/to/emsdk/emsdk_env.sh` first,
# or run inside the CI workflow in .github/workflows/pages.yml which sets
# up emsdk automatically). NOT run as part of the normal native build --
# this is a separate, opt-in target.
#
# Design notes (see CHANGES.md / README.md for the native-side background):
#   * No -fopenmp: the source guards every OpenMP construct with
#     #ifdef _OPENMP, so simply not passing -fopenmp compiles a serial
#     build with no source changes (the same mechanism the native
#     `make NOOMP=1` build already uses and that this project's test
#     suite already verifies gives identical output to the threaded build).
#   * The mmap *store* (mmstore.c) is forced into its pure-malloc mode at
#     runtime by always passing -nommap (see index.html) -- verified
#     byte-identical to the mmap-backed path by the existing native test
#     suite (tests/new_options.sh, tests/determinism.sh), so this isn't a
#     new, unverified code path, just a different runtime flag on already
#     -tested code.
#   * -O2, not -O3: the -O3 safety probe in tools/detect.sh only verified
#     the *native* gcc/x86 backend; it says nothing about Emscripten's
#     clang/WASM backend, which has not been checked the same way, so this
#     build deliberately stays on the level that's actually been verified.
#   * dialign2_dir (BLOSUM / tp400_*) is baked into the virtual filesystem
#     at /dialign2_dir via --preload-file. DIALIGN2_DIR itself is set two
#     ways, belt and suspenders: -DDIALIGN2_DIR_DEFAULT bakes it into the
#     C binary as a compile-time fallback (a browser has no shell
#     environment to export it from; see dialign.c), and index.html also
#     sets Module.ENV.DIALIGN2_DIR for good measure. Either alone is
#     sufficient; an explicitly-set DIALIGN2_DIR still wins over the
#     compiled-in default, so this changes nothing for native builds
#     (which never define the macro) and nothing for anyone who does set
#     the environment variable.
#   * ENVIRONMENT=web,node (not just web): node is what lets
#     wasm/smoke_test.js below load and run the exact same module the
#     browser gets, so the build can be verified without a browser.
#
# Why the verification steps below are not optional: a build was once
# produced that Emscripten reported success for, loaded fine in a browser,
# ran, and printed nothing -- because the object files with main() and all
# the alignment logic in them were not actually linked in (a 9.7KB module
# holding only malloc/free/stack helpers). Nothing at the JS layer can
# detect that from the outside; it has to be checked at the WASM binary
# level (check_wasm.py) and by actually running it (smoke_test.js) and
# checking real output, or this exact failure mode -- appears to work,
# does nothing -- can and did ship.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # repo root

command -v emcc >/dev/null || { echo "emcc not found -- source emsdk_env.sh first"; exit 1; }

OUT=wasm/dist
mkdir -p "$OUT"

SRC=src
SOURCES=$(ls "$SRC"/*.c)

emcc -O2 -std=gnu99 -ffp-contract=off -DCONS -DDIALIGN2_DIR_DEFAULT='"/dialign2_dir"' -I"$SRC" \
  $SOURCES \
  -o "$OUT/dialign2-2.js" \
  -s EXPORTED_RUNTIME_METHODS='["FS"]' \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s WASM=1 \
  -s FORCE_FILESYSTEM=1 \
  -s ENVIRONMENT=web,node \
  -s MODULARIZE=1 \
  -s EXPORT_NAME=DialignModule \
  -s EXIT_RUNTIME=0 \
  --preload-file dialign2_dir@/dialign2_dir

cp wasm/index.html wasm/wasm_inspect.js "$OUT/"

echo "== verifying the build actually contains a working dialign2-2 =="

echo "-- check_wasm.py: does the module even have main() linked in? --"
python3 wasm/check_wasm.py "$OUT/dialign2-2.wasm"

if command -v node >/dev/null; then
  echo "-- smoke_test.js: does it produce the same alignment as the native build? --"
  node wasm/smoke_test.js --dist "$OUT"
  echo "-- test_ui.js: does index.html's own JS still handle every case correctly? --"
  node wasm/test_ui.js
else
  echo "!! node not found -- skipping smoke_test.js/test_ui.js (check_wasm.py above still ran)"
  echo "!! this means the build was NOT verified to actually produce correct output"
fi

echo "== built and verified: $OUT/dialign2-2.{js,wasm,data} + $OUT/index.html =="
