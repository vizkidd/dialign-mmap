# WASM build + GitHub Pages demo

Runs `dialign2-2` entirely client-side in the browser via WebAssembly, at
`https://vizkidd.github.io/dialign-mmap/` once enabled (see "One-time setup"
below).

## What this is
* `wasm/build.sh` -- Emscripten build script. Single-threaded (no `-fopenmp`;
  the source already `#ifdef _OPENMP`-guards every OpenMP construct, so this
  is the same, already-tested code path as the native `make NOOMP=1` build,
  just compiled with a different backend), and always run with the
  equivalent of `-nommap` (the pure-malloc store path, also already verified
  byte-identical to the mmap-backed path by `tests/new_options.sh` /
  `tests/determinism.sh`). `-O2`, not `-O3` -- the `-O3` safety check in
  `tools/detect.sh` only verified native gcc/x86, says nothing about
  Emscripten's clang/WASM backend.
* `wasm/index.html` -- the page: upload a FASTA file (+ optional `.anc`
  anchors), pick common options, run, view/download the `.ali` result.
  Everything happens in-browser; nothing is uploaded to a server.
* `.github/workflows/pages.yml` -- builds with Emscripten and deploys to
  GitHub Pages on every push to `main` that touches `src/`, `wasm/`, or
  `dialign2_dir/`.
* `src/input.c` -- one small, harmless-on-every-other-platform change: the
  input-file `madvise()` call (a pure performance hint, return value
  unused) is skipped under `#ifdef __EMSCRIPTEN__`, since it's not
  guaranteed to exist in Emscripten's libc and there's nothing to gain by
  risking a link failure over a hint. Verified a no-op on the native build
  (full 248-case matrix still identical).

## What is and isn't verified
I don't have an internet connection or Emscripten installed in the
environment I built this in, so **I have not compiled the WASM module
myself.** What I *have* now done, for real, is: (1) execute `index.html`'s actual JS
-- not just read it -- under Node against a real DOM stub
(`wasm/dom_harness.js`) covering every code path (success, silent-no-output,
the new no-`main()` preflight, an `ExitStatus` exit, an unrelated thrown
error, and the anchor-file edge cases), asserted in `wasm/test_ui.js` (20
assertions, all passing, now also run by `wasm/build.sh`), and (2) inspect
byte-for-byte the actual `dialign2-2.wasm` the repo owner built and sent me.

**Update 1:** an actual build hit `callMain is not defined` --
`Module.callMain()` only exists if `EXPORTED_RUNTIME_METHODS` includes it at
build time and evidently wasn't present. Rewritten to use `arguments` +
`preRun` instead -- the same mechanism Emscripten's default auto-run always
uses, nothing extra to export.

**Update 2:** that rewrite then hit `FS is not defined` -- a real scoping
bug: a `preRun` function written in `index.html`'s own script cannot see
`FS`/`ENV` as bare identifiers no matter what's exported, because those are
local to the generated module's own wrapped scope, and closures resolve by
where a function is *defined*, not where it's *called from*. Fixed via the
documented `MODULARIZE=1` behaviour that the config object passed to the
factory *becomes* the Module object in place, so `preRun` now reads
`moduleConfig.FS`/`.ENV` off the object it already holds.

**Update 3 -- the actual root cause of "No .ali output was produced":** I
inspected the `dialign2-2.wasm` from the zip the repo owner sent
(`dialign-mmap-6.zip`, `wasm/dist/dialign2-2.wasm`) at the binary level.
It is 9,713 bytes and defines **10 functions total** -- `malloc`, `free`,
`stackSave`/`stackRestore`/`stackAlloc`, and a few more allocator internals.
**No `main`, no `__main_argc_argv`, nothing from `dialign.c` or any other
source file.** Whatever produced that specific `.wasm` did not actually
link the DIALIGN C sources into it -- the module is valid WASM, loads fine,
`preRun` runs fine (the virtual filesystem exists), and then there is
simply no program to run, so it "succeeds" having done nothing. That fully
explains the symptom: no crash, no thrown error, no log output, just
silence and "done". I can't tell you *why* that particular build produced a
sources-less module from here (stale `dist/` from an earlier placeholder
build? a shell quoting issue that silently dropped `$SOURCES`? worth
checking the raw `emcc` invocation and its exit code next time) -- only
that it did, conclusively, by direct inspection of the bytes.

This is exactly the failure class `wasm/build.sh` now refuses to let
through: it runs `wasm/check_wasm.py` (fails loudly if `main` isn't
exported) and `wasm/smoke_test.js` (fails loudly if the WASM output doesn't
match the native build's, under Node -- see `wasm/golden/`) as part of the
build itself, and the CI workflow inherits that. And `index.html` itself
now does the same `check_wasm.py` check client-side, always (not just in
Debug mode), before ever trying to instantiate the module -- so if this
happens again, the page says exactly this instead of "done, no output".

**Update 4 -- confirmed working, found a real second bug:** a real rebuild
passed `check_wasm.py` (149 functions, `main` genuinely linked in -- Update
3's fix worked) and then `smoke_test.js` caught a second real bug
immediately: `dialign: too many diagonals (more than 0)`. Root cause:
`src/frags.c`'s `frag_init()` sized its address-space reservation with
`(size_t)1 << 32` -- a shift equal to the type's own bit width, undefined
behaviour whenever `size_t` is 32 bits (true on `wasm32`, false on native
x86-64, which is why this never showed up natively). Fixed with an explicit
width check instead of relying on how a given compiler happens to fold
that specific undefined expression -- see `CHANGES.md` for the exact fix
and how it was verified (248/248 native matrix unaffected; the 32-bit code
path itself was force-compiled and real-alignment-tested on this machine,
not just read). This is the pipeline from Update 3 working exactly as
intended: a specific, fixable error instead of silence.

**Update 5:** fixing that revealed the next layer, same failure family:
`dialign: cannot map memory for fragments` on the very first fragment
allocation. The arena reserves address space now (that part works), but
backing part of it with real memory via a second, fixed-address mmap call
does not -- there's no real virtual memory under WASM for that two-step
pattern to work against at all. This one's architectural, not a one-line
fix: `frag_init()` now has a genuine platform split, with the native path
untouched and a WASM-appropriate path that does one plain allocation sized
to one chunk (~63.6 MiB, ~1.85M fragments) instead -- see `CHANGES.md` for
the exact change and how it was verified, including against the precise
`trna100 -n -anc` case from the bug report. If a job needs more than ~1.85M
fragments it will now cleanly hit "too many diagonals" rather than this
error -- that ceiling itself wasn't pushed past in testing (would need a
very large/dense dataset), so it's worth knowing about if a future run
reports that cleanly rather than running to completion.

**Update 6:** confirmed working end to end -- a real run got all the way
through to `anchor.c`'s validation and correctly rejected a mismatched
anchor file (an anchor file from one dataset used with a different
dataset's FASTA file -- check that the two files you pick actually came
from the same source). Found alongside that: a C `exit(1)` does NOT throw
a catchable `ExitStatus` in this Emscripten configuration
(`EXIT_RUNTIME=0`) the way Emscripten's own docs say it does -- the module
promise just resolves normally, same as a clean run, and the only
difference is the missing output file. Nothing was being silently lost
(the full log was already shown), but the page's message undersold it as
generic "no output" rather than naming the likely cause. Fixed -- see
`CHANGES.md` for specifics.

Rebuild with the current `wasm/build.sh` and this specific failure mode
should either not recur, or fail the build itself with a clear message
instead of shipping silently. Still worth double-checking once you rebuild:
1. Does `-anc` correctly find `/input.anc` in the virtual filesystem?
2. Multiple runs on one page load are deliberately *not* supported (each
   click loads a fresh module instance -- see the comment in
   `index.html`) precisely because I could not verify every global the C
   program relies on is safely reset between two runs on the same
   instance. If a fresh-instance-per-run turns out to be too slow in
   practice, that reentrancy would need an actual audit before removing it.
3. Three real bugs have been found and fixed this way so far (`callMain`
   not exported, `FS`/`ENV` not in scope, and now the sources-less build)
   purely from the repo owner pasting the exact browser error or build
   artifact back -- if you hit anything else, that's still the fastest way
   to get it fixed: the actual error or file, not a description of it.

## Debug mode
`index.html` has a "Debug mode" checkbox. When on, every step (file sizes
read, the exact CLI args used, the WASM module's own export table, virtual-
filesystem contents right after the run, full JS stack traces) is logged to
three places at once: the on-page log box, the browser console
(`console.log`/`warn`/`error` as appropriate), and -- since the C program's
own stdout/stderr *are* the `print`/`printErr` callbacks in this
architecture -- what "the module's stdout" means here. The WASM-has-no-`main`
preflight check and the final error message always run/show regardless of
this toggle; Debug mode adds the step-by-step trace on top.

## Verifying a new build yourself
* `python3 wasm/check_wasm.py wasm/dist/dialign2-2.wasm` -- fails if `main`
  isn't linked in. Both `wasm/build.sh` and `index.html` (client-side, via
  `wasm/wasm_inspect.js`, its JS twin) already run the equivalent of this
  automatically; the standalone script is there for manually checking a
  `.wasm` file directly.
* `node wasm/smoke_test.js` -- runs the built module under Node (needs
  `ENVIRONMENT=web,node` at build time, already set) and requires its
  output for two small real datasets to be byte-identical to the native
  build's, stored in `wasm/golden/`. This is what would have caught the
  sources-less build immediately, at build time, with a specific message,
  instead of it reaching a browser at all.

## One-time setup (I can't do this part -- no repo write access from where I ran this)
1. Push this branch/these files to `vizkidd/dialign-mmap` on GitHub.
2. Repo Settings -> Pages -> Source: **GitHub Actions** (one-time toggle).
3. Push to `main` (or run the workflow manually via Actions -> "Build and
   deploy WASM demo to GitHub Pages" -> Run workflow) and watch it build.
4. Site appears at `https://vizkidd.github.io/dialign-mmap/`.
