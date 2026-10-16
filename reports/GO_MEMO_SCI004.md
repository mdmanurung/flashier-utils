# GO memo: what SCI-004a/b still require (flashier-utils)

Date: 2026-10-16. Scope: the two items SCI-004a/b gate a release on, and the
exact evidence needed to call them done. No release without GO on both.

## SCI-004a — verification oracles for the delivered defect fixes

Status at memo time: FND-001 (fresh-native projection oracles) and FND-002
(behavioral partial-fix oracle) are implemented and green in `tests/` —
66 assertions covering both orientations/axes, bit-for-bit reproduction of
from-scratch flashier pipelines (activity/effects projections, PVE, seeded
sampler draws, manifest resets, source_coordinates preservation), and
partial-fixing pinning named entries through full backfit. The multistart
`init_seeds` input-order invariance regression (Oct 8 raw-helper counterexample)
is in place. Full suite: 1230 passed / 0 failures; R CMD check on the built
tarball: 0 errors / 0 warnings / 1 NOTE (optional Suggests absent from the
sandbox, covered in CI).

Remaining for GO:
1. CI green on `main` at `7c214d8` across all engines (ubuntu release/devel,
   macos, windows, full-engines) — observed GREEN 2026-10-16 via check-run
   watch; keep this true at the release commit.
2. A release commit that bumps DESCRIPTION to the release version and records
   the verification receipt (already drafted at the top of HANDOFF.md) in the
   NEWS entry.

## SCI-004b — raw-helper refusal behavior

The Oct 8 seed-consensus input-order counterexample is reconstructed as a
regression test (input-order invariance for `init_seeds`). What remains for GO
is a decision, not code: either (a) the helpers refuse ambiguous input orders
with an explicit error (current behavior if implemented), or (b) they document
the accepted canonical order. GO requires: the chosen behavior stated in the
helper roxygen pages, one test per helper pinning it, and a HANDOFF receipt
line. If neither behavior is yet implemented for any helper, that helper is a
GO blocker.

## Non-goals

No new features beyond the above gate the release; GB-005 (activity threshold +
MDS ordering) is delivered and tested. GBCD/STARCAT remain tracked separately
(see the backlog memo) and do not gate this release.
