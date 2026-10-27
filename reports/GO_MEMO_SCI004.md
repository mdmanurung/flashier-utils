# GO memo: release verification gates REL-001a/b (flashier-utils)

Reconciled 2026-10-27: the two items below were first labelled SCI-004a/b.
That label collides with the frozen protocol
(validation/gbcd-starcat-science/SP-v1/protocol-v1/protocol.md), where SCI-004a
is the fixed-B point projection GO gate (CL-PROJECT-POINT) and SCI-004b is the
supplied-profile consensus recovery GO gate (CL-CONSENSUS). Those scientific
gates remain open and are not addressed here. The software gates in this memo
are renamed REL-001a (verification oracles) and REL-001b (raw-helper input
order); the protocol keeps its IDs.

Date: 2026-10-16. Scope: the two items REL-001a/b gate a release on, and the
exact evidence needed to call them done. No release without GO on both.

## REL-001a (formerly SCI-004a) — verification oracles for the delivered defect fixes

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

## REL-001b (formerly SCI-004b) — raw-helper refusal behavior

The Oct 8 seed-consensus input-order counterexample is reconstructed as a
regression test (input-order invariance for `init_seeds`). What remains for GO
is a decision, not code: either (a) the helpers refuse ambiguous input orders
with an explicit error (current behavior if implemented), or (b) they document
the accepted canonical order. GO requires: the chosen behavior stated in the
helper roxygen pages, one test per helper pinning it, and a HANDOFF receipt
line. If neither behavior is yet implemented for any helper, that helper is a
GO blocker.

Status 2026-10-27: `fit_multistart_ebmf()` is not exported (tests/testthat/
test-core-api.R pins the public API), so no public raw helper exists. Its
`init_seeds` argument documents the accepted order (unique integers, one per
start; R/fit.R) and the input-order invariance regression lives in
tests/testthat/test-fnd-oracles.R. Nothing is exported that needs a refusal
rule. A GO on this item needs no code unless the helper is exported later.

## Non-goals

No new features beyond the above gate the release; GB-005a (activity threshold) and
GB-005b (MDS ordering with canonical sign, declared group order and a dense-
distance refusal) are delivered and tested as of 2026-10-27. GBCD/STARCAT remain tracked separately
(see the backlog memo) and do not gate this release.
