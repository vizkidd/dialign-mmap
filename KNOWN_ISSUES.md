# Known issues

## 1. FIXED: `-fa` output could silently destroy the input file

**Status: fixed in `src/dialign.c` (`would_clobber_input`, `avoid_input_collision`,
`refuse_if_would_clobber_input`), verified below.**

If the input file's name ends in `.fa`, `.seq`, or `.fasta`, DIALIGN strips
that extension internally to build an output-file stem, then reattaches
`.fa` for the `-fa` (separate FASTA-format alignment) output file. When the
input's extension was already `.fa` -- the ordinary, expected naming
convention -- stripping and reattaching reconstructs **the exact same
filename as the input**, and `-fa` output was opened with `fopen(name,"w")`,
which truncates immediately: `./dialign2-2 -anc -n -fa myseqs.fa` silently
overwrote `myseqs.fa` with the gapped alignment output. No warning, exit
code 0. **Confirmed on both the plain original 2.2.1 and this port** before
the fix (identical outcome on both).

### The fix
Before opening any of the four possible output files (`.ali`, `-fa`'s `.fa`,
`-msf`'s `.ms`, `-cw`'s `.cw`), the port now compares the resolved
(`st_dev`,`st_ino`) of the output path against the input file:
* For `-fa` specifically, where the collision is the *common* case (any
  `*.fa`-named input), it does not abort: it writes the alignment to
  `<input>.dialign-aligned.fa` instead and prints a one-line note to
  `stderr`. This keeps `-fa` usable for ordinary filenames while making
  data loss impossible.
* For `.ali`/`.ms`/`.cw`, where a collision should never occur in normal use
  (their suffixes aren't in DIALIGN's extension-stripping list), it refuses
  with a clear error rather than silently truncating, as a defensive
  backstop.

This is a deliberate, intentional divergence from the original's behaviour
for exactly the input-destroying case, and only that case -- every other
output byte for byte, for every one of the 248 cases in
`tests/cases.txt`, is unchanged (`tests/compare.sh` was updated to
recognise and pass this specific, expected difference explicitly, labeled
`ORIG-CLOBBERS-INPUT(fixed);NEW-PRESERVES-INPUT`, rather than silently
special-casing it).

### Verification
* `tests/run_matrix.sh` against a clean original build: **ALL 248 CASES
  IDENTICAL** (3 of these are the fixed `-fa` collision cases, now passing
  under the new, explicit `ORIG-CLOBBERS-INPUT` label instead of failing as
  a spurious byte-diff; another 3 are the pre-existing, already-documented
  `-ma`/`-lgsx` `ORIG-CRASHES` cases from CHANGES.md section 5b).
* Manually confirmed: `-anc -n -fa` on a `.fa`-named input no longer
  modifies the input file (byte-identical before/after, 0 gap characters
  before and after), writes the correct alignment to
  `<input>.fa.dialign-aligned.fa`, and exits 0.
* Manually confirmed: `-fa` on a non-colliding filename (e.g. `foo.seqdata`)
  is completely unaffected -- same output filename and content as before
  the fix.
* Full fuzz pass (60 runs, seed 1) after the fix: same, already-diagnosed
  `-ma` mismatch pattern only, nothing new.

## 2. `tests/fuzz.py` runs each binary twice (`r`/`r2`, `n`/`n2`) but only compares `r` to `n`
The second pair of runs is computed and discarded; nothing checks `r==r2` or
`n==n2` (which is presumably what they were for -- a same-process-run
determinism check). Not a bug in `dialign2-2` itself, just dead/incomplete
code in the fuzz harness. Still open; low priority.
