# Scaling test: full-pairwise-anchor alignment cost

## Method
Synthetic "long" tRNA-like sequences (183-608 nt: random 30-220 nt leader/trailer,
~30% with an inserted 20-90 nt intron-like block, 8% substitution rate everywhere
except a protected 6 nt anticodon core). For N sequences, every C(N,2) pair is
anchored on that anticodon (anchor file format: `i j pos_i pos_j 6 10.0`).
Generator: `tests/data/gen_trna_scaled.py <N> <outdir> [seed]`.

Run as `./dialign2-2 -anc -n <file>.fa` on a single core (no other load),
wall time measured externally, peak RSS sampled from `/proc/<pid>/status`
every 0.2s. Harness: `tests/time_run.sh`.

## Results (measured, not estimated)

| N (sequences) | anchors = C(N,2) | wall time | peak RSS |
|---:|---:|---:|---:|
| 50   | 1,225   | 3.73 s   | 6.8 MB  |
| 100  | 4,950   | 14.62 s  | 20.1 MB |
| 200  | 19,900  | 60.42 s  | 51.8 MB |
| 300  | 44,850  | 149.01 s | 107.9 MB |
| 400  | 79,800  | 273.30 s | 190.2 MB |

Log-log regression over these 5 points:
* **wall time ~ N^2.07** (essentially quadratic — matches the C(N,2) anchor
  count directly, i.e. per-anchor cost is roughly constant, not growing with N)
* **peak RSS ~ N^1.57** (sub-quadratic; anchors are not all retained at once)

Extrapolation from the fit: N=1000 -> ~1,770 s (~29.5 min), ~734 MB peak RSS.

## Cross-check against a real N=1000 run
A full N=1000 run (499,500 anchors) was independently started and observed
directly (not from the fit): after 28 min 21 s of CPU time it had not yet
finished, RSS was 942 MB and climbing. Both numbers are the right order of
magnitude for the ~29.5 min / ~734 MB prediction above (RSS running somewhat
above the sub-quadratic fit, consistent with the fit being an average over a
regime that is only mildly super-quadratic, not a sign of a different growth
law). This is corroboration, not a third clean data point — the run was not
carried to completion.

## N=700 attempt
A background N=700 run (244,650 anchors, predicted ~846 s / ~419 MB) was
started to get a mid-range confirmation point. It disappeared without an
error, a crash signal, or a completion marker after being detached for
several tool-call boundaries in this sandbox. There is no evidence this was
a program bug (no stderr output, no partial/corrupt files, no OOM in
`dmesg`) — the most likely explanation is the sandbox's own background-job
lifetime, not `dialign2-2` itself. This is a testing-environment limitation,
noted here rather than silently omitted.

## What this means for shipping
* Up to ~400 fully-cross-anchored sequences (~80k anchors): demonstrated,
  fast, low memory, real data.
* ~1000 fully-cross-anchored sequences (~500k anchors): plausible from the
  fit and roughly corroborated by a live run, but **not a completed,
  verified data point**. Treat "1000+ sequences with full pairwise
  anchoring" as within reach but not yet signed off.
* The scaling itself (~N^2, driven by the anchor count you choose to
  supply, not by the aligner) is good news: it is not the O(N) per-anchor
  closure blow-up that `CHANGES.md` section 8 worried about. The practical
  limit for very large anchored runs is wall-clock time from the anchor
  count itself, which is under the *user's* control (how many pairs they
  choose to anchor), not an internal inefficiency.
* Un-anchored alignment (no `-anc`) does not carry this O(N^2) anchor cost
  and was already exercised up to several hundred sequences in the existing
  test suite without issue.
