# dialign-mmap: changes, fixes, improvements

## 0. Why the code base was re-based

The previous dialign-mmap kept its working data in *text* files (`.diags`,
`.op`, `.ow`, `.ancwork`) addressed by line number and cached byte offset.
Tested against a build of the original DIALIGN 2.2.1 it was not correct:

* `test1.fa` (3 seqs): different name column width (`SEQ_NAME_LEN` 50 vs 12), different
  weight row, tree with `inf` branch lengths; `test2.fa` different and 8x slower;
  `trna_50.fa` segfault after 195 s (original: 5 s).
* weights/sums had been changed from `float` to `double`; the original's non-stable
  quicksort and tie handling depend on the exact float values.
* working files were left behind next to the input.
* known-open bugs: >35 sequences, `-anc`, an `-O3` miscompile worked around with `-O0`.

The new tree starts from the original sources (so algorithms, float types and
tie-breaking are the original's) and moves the memory into an mmap store.
The `-O0` workaround, the `-anc` name collision and the >35-sequence failure
disappear by construction.

## 1. The mmap store (`mmstore.c/h`)
* size-class allocator over 64 MiB file-backed chunks (per-thread free lists, no lock on the
  fast path); blocks >= 1 MiB get their own file mapping and are `munmap`ed on free.
* files are anonymous (`O_TMPFILE`, or `mkstemp`+`unlink`): nothing is left behind, even after a crash.
* fresh pages read as zero: `mm_calloc` costs nothing; untouched pages use no memory or disk.
* `SIGBUS` handler explains a full temp volume; if temp files cannot be created it falls
  back to the heap with one warning.

## 2. What now lives in the store
sequences and names (no MAX_SEQNUM, no line-length limit); `amino`/`amino_c` (now 1 byte per residue instead of 4);
`open_pos` (one contiguous block, 1 byte per flag instead of 4 and stored *inverted* so that the
N^2*L table needs no initialisation and untouched parts use no memory/disk);
`glob_sim`, `cont_it_p`; the whole alignment-graph closure (`pratique.c` allocators);
per-pair DP arrays and candidate fragments; all fragment lists and sort arrays;
output buffers of `ali_arrange` (`endseq`, `hseq`, `inv_shift`, `shift`, fragment copy);
`-mat` substitution table (now one flat block instead of N^2*21 allocations); `-xfr` exclusion lists; motif positions;
UPGMA similarity matrix.
Dropped: `pair_score` (N x N floats that the original writes but never reads).

## 3. Speed
* `print_log()` (`-lo`/`-fsm`/`-fsmv`) walked the entire fragment list once per sequence pair (O(N^2 * D)), and was
  called even when it printed nothing. It is now skipped when nothing is requested, and otherwise the diagonals are
  bucketed per pair in one stable pass (O(N^2 + D)); the `FRG n` numbers of the original counting rule are reproduced.
  200 tRNAs: 69 s -> 13.6 s (no log options).
* `av_tree_print`: cached best partner per row (same tie-breaking), dynamic clade names
  (the original reserved `100*N` bytes for every clade = O(N^2) memory), matrix in the store.
* pairwise stage parallelised with OpenMP (see 4).

## 4. Parallel pairwise alignment (deterministic)
Pairs are handled in blocks of 512; inside a block they run in parallel and the results are appended to the
fragment list in the original (i,j) order. Races removed: the global weight tables that `rel_wgt_calc`
rewrites per pair (now per-thread tables; `rel_wgt_calc` takes an explicit type instead of comparing
table pointers), the global `pair_dia` (returned per call), the counters `dia_num`, `max_dia_num`, `cont_it`
(atomic). `ow_add()` reads the weight tables left behind by the *last* pair, so they are re-created
after the parallel stage. Serial fallback for options that write files from inside the stage.

## 5. Bug fixes (original behaviour that was wrong or undefined)
* `throw_out()` freed the first rejected fragment of an iteration without re-linking its predecessor
  (dangling pointer). Fixed.
* `word_count()`: `strlen(s)-1` underflow on an empty string; last character of a line without newline ignored
  (breaks `.anc` files without final newline). Fixed.
* FASTA reader: heap over-run by one byte in `seq_shift`, wrong buffer sizes for files without final newline,
  `\r` kept in names for CRLF files, 10 000-character line limit. Fixed (well-formed input: unchanged behaviour).
* `-mat` together with `-n` crashed (segfault in the original); now a clear error message.
* `erreur()` called with 2 arguments; missing prototypes; `int main`; makefile: `-ffp-contract=off`
  (bit-identical floats on FMA platforms), `-std=gnu99`.

## 5b. Further bugs found by the option sweep / robustness tests (all present in the original)
* **Uninitialised variable in `frag_chain`** (`start_pep_c` was read even without `-cs`; the value depended on stack
  garbage, so results changed with unrelated things such as `-pst`). Found with `-ftrivial-auto-var-init`
  builds; initialised, and the unguarded read is now guarded. `-lgs_t`, `-it ... -pst` now match the original.
* **Crash with mixed alignments** (`-ma`, `-lgsx`, `-nt -ma`): in mixed mode a translated fragment can end up to
  3 positions behind the sequence end, and the original overran `diap[]`/`prec_vec[]` (segfault on e.g. 20 DNA
  sequences). The arrays now have `PADX=8` slack; inputs that did not crash give unchanged results (the padded
  slots are never read by the sweep — `for(i=1;i<=seqlen[n1];i++)` only ever consumes `diap[i+1]`, so a
  candidate landing past `seqlen[n1]+1` was always dead weight, not a result-affecting one). The post-sweep
  free loop was widened from `hv<seqlen[n1]+3` to `hv<seqlen[n1]+PADX` to match, so the one slot a translated
  candidate can actually land on (`seqlen[n1]+3`) doesn't leak a `pair_frag` node. Re-verified against a plain
  (non-ASan) rebuild of the unmodified original: 248/248 matrix cases identical, including the 3 cases where
  the original segfaults and the port does not; 200/200 fuzz runs (seeds 1, 2) with no unexplained mismatches;
  248/248 robustness cases stable across all 7 build/run variants.
* **Closure tables** are zero-filled (the original read never-written heap memory and only worked because fresh heap
  pages happen to be zero; the original's output changes under `MALLOC_PERTURB_=165`).
* `throw_out()` re-links via a correctly maintained tail pointer.
* `DIALIGN_MMDEBUG=1` poisons freed store blocks (use-after-free detector).
* **`-fa` could silently destroy the input file** (present in the original too): when the input's own name ends in
  `.fa`/`.seq`/`.fasta`, the stripped-and-resuffixed `-fa` output filename can reconstruct the exact input filename,
  and `fopen(...,"w")` truncated it on the spot -- `./dialign2-2 -anc -n -fa myseqs.fa` silently replaced
  `myseqs.fa`'s raw sequences with the gapped alignment, no warning, rc=0. Fixed with an (`st_dev`,`st_ino`)-based
  collision check: `-fa` now falls back to `<input>.dialign-aligned.fa` with a one-line stderr note (the common,
  `*.fa`-named-input case); `.ali`/`.ms`/`.cw` hard-refuse instead (should never trigger in normal use, defensive
  only). This is a deliberate, intentional divergence from the original for exactly this input-destroying case only
  -- re-verified against a plain rebuild of the unmodified original: 248/248 matrix cases identical (3 of these are
  the fixed `-fa` cases, explicitly labeled `ORIG-CLOBBERS-INPUT` in `tests/compare.sh` rather than silently
  special-cased or left as a spurious byte-diff), fuzz pass unaffected. See `KNOWN_ISSUES.md` #1.
* **Uninitialized `fp_ali` crash/hang with `-pand`/`-o` + a textual-alignment-disabling mode** (`-lgs`, `-lgs_t`, or
  `-nta`): `fp_ali` is only `fopen()`'d when `textual_alignment` is set, but the `-pand`/`-o` (`pr_av_nd`/
  `pr_av_max_nd`) summary lines were `fprintf(fp_ali, ...)` unconditionally, so e.g. `-lgs_t -pand` wrote through
  an uninitialized stack `FILE*` -- undefined behaviour, reproduced as a deterministic SEGV under ASan and, on
  an unmodified build, as a *non-deterministic* crash/timeout across repeated runs of the identical input
  (garbage stack contents vary run to run) -- exactly what a fuzz run surfaced (`-nt -cw -o -max_link -pand
  -lgs_t`: original timed out on one run and segfaulted on an identical repeat). Fixed by gating both
  `fprintf(fp_ali, ...)` calls on `textual_alignment`, matching the `fclose(fp_ali)` right after them which
  already did. No effect when `textual_alignment` is set (the ordinary case). Confirmed present in the original
  2.2.1 too. Re-verified: 248/248 matrix identical; fuzz seeds 1-2 (200 runs) show only the already-catalogued
  `-ma`/`-lgsx` original-crashes pattern, including this exact case and one previously-undetected sibling
  (`-pamnd -nta`, same root cause, same fix) -- no more timeouts or unexplained mismatches anywhere.

## 5c. More parallelism
* Overlap weights (`ow_add` for all diagonal pairs, O(D^2)) are computed in parallel: every diagonal accumulates its
  own contributions in the order the serial loop would have added them, so the float sums (and output) are
  bit-identical for any thread count; serial code path kept for 1 thread and for `-pst`.

## 5d. Fragment memory, sequence-count limit, guide tree
* **Compact fragments.** `struct multi_frag` (one per selected/candidate diagonal) shrank from 56 to 36 bytes:
  the two 8-byte list pointers became one 32-bit index (`next`; `pred` was only used by a workaround and is gone),
  `it` is 16 bit, `sel/trans/cs` are bit-fields. The nodes live in one file-backed arena (`frags.c`: a reserved
  address range that is filled chunk by chunk with temp-file mappings, so an index maps to an address by one
  multiplication). Limits: 2^32-1 fragments, at most 65 535 iteration steps (`-istep`). The arena is only touched
  from serial code (list building, `throw_out`, anchors); the parallel regions never allocate list nodes.
* **No MAX_SEQNUM.** The original has `#define MAX_SEQNUM 10000` and static arrays (built with 200 it segfaults on 400
  sequences, and it segfaults on 10 050). All per-sequence arrays are dynamic now. Checked: byte-identical to the
  original built with a larger MAX_SEQNUM on 400 sequences; the reader accepts 10 500 sequences
  (`tests/new_options.sh`). A *complete* alignment of > 10 000 sequences was NOT demonstrated: the closure update
  costs O(N) per accepted diagonal (as in the original), which is hours on the single core available for testing.
* **Guide tree.** `av_tree_print`: row-maximum cache (see 3); rows are not rescanned when the best partner's
  similarity did not decrease, which removes an O(N^3) behaviour for all-equal similarities. Average linkage is still
  O(N^3) in the worst case (as in the original). New option `-notree` skips the tree (the tree line of the output stays empty).
* `-ma` / `-lgsx`: see 5b (padding). The mixed-mode design flaw (translated fragment candidates are recorded at
  codon-unaligned lengths) is deliberately not changed, because that would change results of runs that did not crash.

## 5e. Top-level build system (`Makefile`, `tools/detect.sh`)
A root `Makefile` wraps `src/makefile` with `build`/`install`/`uninstall`/`test`/`test-full`
targets. `make` auto-detects OpenMP support (a real compile+link probe, not just checking
the compiler name) and whether `-O3` is safe on this compiler/codebase: it builds both
`-O2` and `-O3`, runs both on a few real datasets, and only uses `-O3` if the output is
byte-for-byte identical (the `-ffp-contract=off` guarantee this whole project depends on
means "the flag is accepted" is not the same as "safe to use" -- see `tools/detect.sh`).
Verified: full 248-case matrix identical between an `-O3` build and the plain original,
not just `-O3`-vs-`-O2` self-consistency. `make install` places a small wrapper script
(sets `DIALIGN2_DIR`, execs the real binary) alongside the reference data, so the tool
works immediately after install with no manual environment setup -- the single most
common first-run error with the plain original.

## 5f. WASM / GitHub Pages build (`wasm/`, `.github/workflows/pages.yml`)
An Emscripten build (`wasm/build.sh`) targeting the browser: no `-fopenmp` (the source already
`#ifdef _OPENMP`-guards every OpenMP construct, so this is the same code path as the native `NOOMP=1`
build, just a different compiler backend) and always run with the equivalent of `-nommap` (already
verified byte-identical to the mmap-backed path natively). One source change: `src/input.c`'s
`madvise()` on the mmap'd input file (a pure hint, return value unused) is skipped under
`#ifdef __EMSCRIPTEN__`, since it isn't guaranteed to exist there and there's nothing to gain by
risking a link failure over a hint -- confirmed a no-op on native (248/248 matrix unaffected).
**Not yet compiled or run in an actual browser** -- built without internet access or Emscripten
available; see `wasm/README.md` for exactly what is and isn't verified and what a first real build
needs to check. **Update 1:** the repo owner did build it and hit `Module.callMain is not defined`
(it needs to be added to `EXPORTED_RUNTIME_METHODS` at build time and evidently wasn't present in
the built `dialign2-2.js`). Rewritten to use Emscripten's `arguments`+`preRun` config keys instead
of `Module.callMain()` -- the same mechanism the default auto-run startup path always uses, so
nothing extra needs exporting for it to exist. **Update 2:** that then hit `FS is not defined` --
an actual scoping bug: a `preRun` function defined in the page's own script cannot see `FS`/`ENV`
as bare identifiers regardless of exports, since those are local to the generated module's own
wrapped scope. Fixed via the documented `MODULARIZE=1` behaviour that the config object passed to
the factory function becomes the Module object in place, so `preRun` now reads
`moduleConfig.FS`/`.ENV` off the object it already holds a reference to.
**Update 3 -- root cause found:** the next report was "done, no .ali output, no log at all" (silent
success). Direct byte-level inspection of the repo owner's actual `dialign2-2.wasm`
(`wasm/check_wasm.py`, `wasm/wasm_inspect.js`) showed 10 defined functions total -- malloc/free/
stack helpers only, no `main`/`__main_argc_argv`, none of the DIALIGN C sources linked in at all.
The module was valid WASM and loaded fine; there was just no program in it to run. `wasm/build.sh`
now runs `check_wasm.py` and `wasm/smoke_test.js` (runs the built module under Node, requires its
output to match the native build's for two golden cases) as part of the build itself and fails
loudly if either doesn't pass; `index.html` runs the same check client-side, always, before
instantiating the module, so this exact failure mode now gets a specific error instead of silence.
Also added: a `DIALIGN2_DIR_DEFAULT` compile-time macro (`dialign.c`, WASM build only, native
builds never define it and are unaffected -- confirmed, 248/248 matrix unchanged) so the module
has a working default even if `Module.ENV` isn't reached; a Debug-mode toggle in `index.html`
(logs args/file sizes/module exports/VFS contents/stack traces to the page, the browser console,
and the module's own stdout/stderr at once); and `wasm/dom_harness.js` + `wasm/test_ui.js`, which
actually execute (not just read) `index.html`'s JS under Node against a stubbed DOM across every
branch (success, silent-no-output, no-main, ExitStatus, unrelated throw, anchor-file edge cases)
-- 20/20 assertions passing, now run by `build.sh` too.
**Update 4 -- real build, real bug found by the new pipeline:** the repo owner rebuilt with a
working toolchain. `check_wasm.py` now passes (`main` genuinely linked in, 149 functions) --
confirms Update 3's fix worked -- but `smoke_test.js` caught a second, real bug immediately:
`dialign: too many diagonals (more than 0)`. Root cause: `frags.c`'s `frag_init()` sized its
upfront address-space reservation with `size_t slots = (size_t)1 << 32`, a shift equal to the
type's own bit width -- undefined behaviour whenever `size_t` is 32 bits. Native x86-64 builds
never hit this (`size_t` is 64 bits there); `wasm32` has a 32-bit `size_t` and evidently folds
this to a value that makes every reservation attempt fail immediately, silently, with no
address-space-exhaustion path taken at all. Fixed with an explicit `#if SIZE_MAX > 0xFFFFFFFFu`
split: the exact original target (2^32 nodes, same floor) on any platform wide enough to express
it without UB, and a smaller, still-generous, explicitly-chosen opening bid/floor (16M/64K nodes)
on narrower ones, instead of leaving the result to compiler-specific UB-folding. Verified: 248/248
native matrix unaffected; the `#else` branch was force-compiled on this machine (`-USIZE_MAX
-DSIZE_MAX=0xFFFFFFFFu`, real build, not a simulation of 32-bit integer semantics, but of this
specific code path) and produces byte-identical output to the normal build on real alignments,
including one (`trnaL_80 -anc`) large enough to actually exercise the second-chunk-allocation path
this function guards. The pipeline from Update 3 is working as designed: a real, specific,
diagnosable error from an automated check, not a silent failure reaching a browser.
**Update 5 -- third real bug from the same pipeline, same root cause family:** fixing Update 4
revealed the next layer: `mm_reserve()`'s `PROT_NONE` reservation now succeeds on WASM (Emscripten's
mmap shim appears to treat it as a no-op, consistent with reservations being "free" by design), but
the *second* step -- `mm_map_fixed()`'s `MAP_FIXED` mmap to back part of that same reservation with
real memory at that exact address -- does not: `dialign: cannot map memory for fragments` on the
very first fragment allocation. That reserve-now/commit-later-at-a-fixed-address pattern has no real
equivalent in a single growable linear-memory model; there is no actual OS-level virtual memory
under WASM for a second mmap call to "fix up" an earlier reservation. This is architectural, not a
typo, so the fix is a genuine platform split in `frag_init()` rather than a one-line correction: the
`SIZE_MAX > 0xFFFFFFFFu` (native) branch is untouched, byte-for-byte; the narrower-`size_t` branch
now skips the reserve/commit scheme entirely and does one ordinary `calloc()` sized to exactly one
chunk (`CHUNK_BYTES`, ~63.6 MiB, ~1.85M fragments -- the same per-chunk size the native path grows
by), with `mapped` set equal to `max_slots` from the start so a job needing more than that takes the
existing, already-tested "too many diagonals" exit instead of ever reaching the broken
`mm_map_fixed()` call again. Verified: 248/248 native matrix unaffected (that branch is literally
unchanged code); the `#else` branch force-compiled and real-alignment-tested on this machine again,
byte-identical to the normal build on `prot6`, `dna8 -n`, and -- matching the exact scenario from
the bug report -- `trna100 -n -anc`. The fragment-count ceiling itself (1.85M) was not empirically
pushed past in this round (would need a pathologically large/dense dataset to construct); confidence
there rests on the arithmetic (`mapped == max_slots` makes the overflow check fire immediately,
routing into exit logic already covered by existing tests) rather than a new empirical test, and
that's flagged here rather than left implicit.
**Update 6 -- Update 5 confirmed working; a UX/correctness-of-messaging gap found alongside it:**
a real run (`trna_100.fa` paired with `trna100.anc` -- two different datasets that happen to have
similar names, not the same dataset's own anchor file) got all the way through to `anchor.c`'s
anchor-vs-sequence-length validation and correctly rejected it with `PROGRAM TERMINATED`, confirming
Update 5's fix works end to end. But the page showed the generic "No .ali output was produced"
message rather than the more specific "exited with code N" one, because of something worth
recording: `index.html`'s comments claimed a C `exit(1)` "always" surfaces as a thrown `ExitStatus`
in Emscripten -- true per Emscripten's own docs, but a real build with `EXIT_RUNTIME=0` was observed
to NOT throw there, resolving the module promise exactly as if `main()` had returned normally. The
only visible difference is the missing output file, which the code already checked for and already
showed the full log for -- so nothing was silently lost, but the messaging undersold it as a generic
"no output" rather than flagging it as the likely-error case it almost always is. Fixed: the
no-output message now explicitly names the common cause (anchor file not matching its FASTA file)
and treats a missing output file as the primary failure signal rather than a fallback for when the
`ExitStatus` catch "doesn't happen to fire" -- it's the common case, not the edge case, for this
build configuration. `wasm/test_ui.js`'s `noop` scenario (and its assertions) updated to match this
confirmed-real behaviour; the `ExitStatus`-throwing scenario is kept as a defensive secondary case
since it may still apply to some Emscripten version/config, just not the one actually observed.

## 6. Deliberately unchanged original quirks
`-red`/`-strict` are rejected by the argument check; `ow_add` uses the last pair's weight table;
`SEQ_NAME_LEN` is 12 (override with `-DSEQ_NAME_LEN=n` if you want longer names, at the cost of a different layout).

## 7. Verification (all against a build of the original 2.2.1 source, byte comparison of every output file)
* `tests/run_matrix.sh`: 250 dataset x option cases from `tests/cases.txt` (protein, DNA, `-n -nt -ma -lgs -cs -ds -o -ow -iw -it -anc -xfr -sc
  -ref_seq -mot -mat -stdo -fn -fa -msf -cw -ff -fop -afc -afc_v -lo -fsm -bs -pao -csc -mask -max_link -min_link ...`,
  2..50 sequences, tie-heavy sets, messy FASTA, 5 kb sequences): ALL IDENTICAL.
* 100, 200 and 300 tRNA sequences (`-n`): identical (300: original 475 s, this version 28 s).
* `tests/determinism.sh`: 1/2/3/4/8 threads, `-nommap`, `-tmpdir`, no leftover files: OK.
* AddressSanitizer + UBSan build over the main paths: no reports.
* Resident anonymous memory on 100 tRNAs: original 7.3 MB, mmap build 0.2 MB.

* Anchored alignments (`-anc`), all identical to the original: 100 short tRNAs (4950 anchors, also `-nas`, protein mode),
  300 short tRNAs (44 850 anchors; original 665 s, this version 38 s), 40 and 80 *long* tRNAs (300-600 nt with random
  leader/trailer and optional introns; 780 / 3160 anchors). Every anchored residue pair is verified to sit in the same
  alignment column (29 700 of 29 700 for the 100 short set, 4 680 of 4 680 for the 40 long set).
* `tests/robustness.sh`: every case must give the same output for default / `-nommap` + `MALLOC_PERTURB_` / poisoned frees /
  `-ftrivial-auto-var-init=zero` and `=pattern` / `-O0` / `-threads 3`.

* `tests/fuzz.py`: random inputs x random option sets against the reference (original + minimal bug fixes); also checks that
  each program is repeatable. `tests/tree_unit.sh`: original vs new guide-tree code on 1 200 random tie-heavy matrices.
  `tests/new_options.sh`: `-nommap -threads -tmpdir -notree`, unusable tmpdir, bad arguments, 10 500-sequence reader.

## 8. Not done / limits (please read)
* Only a 1-core machine was available: thread-count *determinism* was tested (up to 8 threads time-sliced),
  **speed-up was not measured**; ThreadSanitizer was not run.
* Overlap weights are still O(D^2) work (now parallel, see 5c).
* Fragment nodes are 36 bytes (was 56): about 36 % less, not half.
* Cost of the closure update per accepted diagonal is O(N) (original algorithm) - this, not memory, limits very large N.
* A completed alignment of more than 10 000 sequences has not been demonstrated (see 5d).
* Two mismatches seen in the very first fuzz run (seed 1, iterations 27 and 119) could not be reproduced and are unexplained;
  more than 1 000 later fuzz iterations with several seeds had none.
