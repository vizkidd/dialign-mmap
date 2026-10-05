# Recommended next features / hardening (not blockers, but worth planning)

Grounded in what actually came up while testing this port, not a generic wishlist.

1. **Progress/verbosity for long anchor-heavy runs.** A `-anc` run on a few
   hundred+ sequences gives zero feedback for minutes at a time — no way to
   tell "working normally" from "stuck." Even a coarse `stderr` line per
   N processed anchors would remove the guesswork this review needed
   background polling and `/proc` sampling to work around.

2. **A pre-flight anchor-count warning.** Anchor cost scales ~N^2 with the
   number of anchors supplied (see PERFORMANCE.md) — that's fine and
   expected, but a user who mechanically anchors every pair for a
   1000+-sequence set has no way to know that means ~30 min and ~500k
   anchor lines until they've already started. A one-line estimate printed
   at start (`~N anchors, expect ~T seconds`) based on the measured
   ~N^2.07 curve would let users choose sparser anchoring (e.g. a spanning
   tree of anchors instead of all pairs) knowingly instead of by accident.

3. **DONE: the `-fa`/input-filename collision (KNOWN_ISSUES.md #1) is fixed** —
   was the single highest-risk item found in this review (silent data loss);
   now auto-avoided for `-fa` and hard-refused for `.ali`/`.ms`/`.cw`, verified
   against the full 248-case matrix and a fuzz pass.

4. **ThreadSanitizer pass.** `CHANGES.md` section 8 already flags this as
   not done; the OpenMP pairwise stage is exactly the kind of code where a
   clean ASan/UBSan run (already done) doesn't rule out data races that TSan
   would catch. Worth doing before calling the parallel path fully verified.

5. **A completed, verified large-N anchored benchmark.** This review got a
   clean, real curve up to N=400 and corroborating-but-incomplete evidence
   at N=700/1000 (PERFORMANCE.md). Actually finishing one clean run in the
   1000-2000 range, on a machine that can hold a background job for the
   full duration, would close the one open scaling question this port
   currently has.

6. **Sparse/streaming anchor format for very large anchor sets.** Right now
   a full pairwise anchor file for N sequences is O(N^2) lines
   (500k lines / 11 MB at N=1000) generated and parsed as plain text. If
   large fully-anchored runs become a real use case, a documented
   alternative (e.g., "anchor only a spanning set of pairs, not all of
   them") would sidestep the scaling question in (5) entirely rather than
   needing to make full pairwise anchoring itself faster.

7. **Wire up or delete the dead `r2`/`n2` determinism check in `fuzz.py`**
   (KNOWN_ISSUES.md #2) — cheap, and determinism-under-repeat-run is
   exactly the kind of thing a fuzzer should be checking.

None of these are correctness regressions in the port itself — the
port's actual alignment logic has held up under everything this review threw
at it (matrix, fuzz, robustness, and the new scaling test). They're gaps
in tooling, documentation, and one inherited (not introduced) footgun.
