# dialign-mmap

DIALIGN 2.2.1 with all size-dependent data structures kept in mmapped
temporary files (low resident memory, no limit on sequence length or count
other than disk space) and an OpenMP-parallel pairwise-alignment stage.
**Output is byte-identical to the original DIALIGN 2.2.1 for the same
command line, for every thread count.**

**Browser demo (WASM):** see `wasm/README.md` -- once enabled, a
single-threaded WASM build runs entirely client-side at
`https://vizkidd.github.io/dialign-mmap/`. Not yet compiled/tested in an
actual browser as of this writing (built in an environment without
Emscripten or internet access) -- see that file for exactly what's
verified vs. what still needs a first real build to confirm.

## Build
    make                       # auto-detects OpenMP and whether -O3 is
                               # safe (byte-identical output, not just
                               # "compiles") on this compiler; falls back
                               # to -O2 otherwise -- see tools/detect.sh
    make CC=clang              # use a specific compiler
    make NOOMP=1               # force serial (skips OpenMP detection)
    make OPT=-O2               # force an optimization level (skips the -O3 probe)
    make print-config          # show what was detected/chosen, without building

Lower-level, same options, no auto-detection or install/test targets:

    cd src && make             # OpenMP + mmap store, -O2
    make DEBUG=1               # -O0 + AddressSanitizer/UBSan

## Install
    make install                          # -> /usr/local
    make install PREFIX=$HOME/.local      # -> a user prefix
    make install DESTDIR=/pkg PREFIX=/usr # staged install, e.g. for packaging

Installs a small wrapper script named `dialign2-2` that sets `DIALIGN2_DIR`
for you and execs the real binary (`dialign2-2.bin`), plus a copy of the
reference data it points at -- so after `make install` you can just run
`dialign2-2 seqs.fa` with no environment variable to set up first.
`make uninstall` removes exactly what `install` put down.

## Run
    export DIALIGN2_DIR=$PWD/dialign2_dir   # only needed if run in place,
    ./src/dialign2-2 [options] seqs.fa      # not after `make install`

All original options are unchanged. New options:

| option          | meaning |
|-----------------|---------|
| `-tmpdir <dir>` | where the mmap store puts its (unlinked) temp files. Default: `$DIALIGN_TMPDIR`, else the directory of the input file. **Use a disk, not tmpfs.** |
| `-nommap`       | keep everything in ordinary heap memory (fast for small inputs, same output) |
| `-notree`       | skip the guide tree (average linkage is O(N^3) worst case; the tree line of the output is left empty) |
| `-threads <n>`  | OpenMP threads (default: `OMP_NUM_THREADS` / all cores). Options that write per-fragment side files from inside the pairwise stage (`-afc -fop -pst -lo -mot -it -wgtpr -pamnd`) run that stage serially, so their files are identical too. |

## Tests
    make test                                      # self-contained: new_options.sh
                                                    # (+ tree_unit.sh if ORIG_SRC given)
    make test ORIG=/path/to/original/dialign2-2 \
              ORIG_SRC=/path/to/original/src        # + full differential suite:
                                                    # matrix, determinism, fuzz

Lower-level, run individually against a build of the original 2.2.1 source
(`dialign-2_2_1-src.tar.gz`, built separately -- it is not part of this repo):

    tests/run_matrix.sh  /path/to/original/dialign2-2  src/dialign2-2   # 248 dataset x option cases
    tests/determinism.sh /path/to/original/dialign2-2  src/dialign2-2   # threads / mmap / tmpdir
    tests/robustness.sh src                                              # UB / uninitialised-memory detector (7 build variants)
    tests/fuzz.py <reference-bin> src/dialign2-2 200 1                   # differential fuzzing
    tests/tree_unit.sh /path/to/original/src/functions.c
    tests/new_options.sh src/dialign2-2
See CHANGES.md for the full list of changes, KNOWN_ISSUES.md and
FEATURES.md for what's still open.
