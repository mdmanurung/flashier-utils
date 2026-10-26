# StarCAT status

Date: 2026-10-25. Scope: this flashier-utils checkout.

The dedicated cross-dataset StarCAT catalog/reference pipeline remains planned.
The revised design bundle exists in flashier-starcat/. Its execution tracker
contains TODO core tasks and deferred experimental prior-transport tasks.
Current source and exports contain no shipped factor_reference catalog builder,
reference consensus graph workflow or completed cross-dataset reference study.

## What exists

- Matched factor representations and extraction, loading congruence/alignment,
  native fixed-factor projection, holdout and stability primitives.
- The exact-pin GBCD adapter and its recorded installed-engine qualification.
- A revised StarCAT design that separates generic comparison profiles from
  source-faithful StarCAT profiles and prevents assignment from pruning graph
  edges. It also separates residual-program discovery from the core catalog.

These foundations do not complete the dedicated StarCAT workflow.
The current cNMF pilots explicitly use build_ref=False; they created no starCAT
reference export. The GBCD comparison's common gene-association profiles are
also a private descriptive endpoint, not a StarCAT implementation.

## What remains

SC-001: source collection/eligibility and validated comparison profiles.
SC-002/003: all-pairs similarities, graph candidates and separate generic versus
StarCAT merge rules. SC-004/005: consensus schema, stable identities, reference
builder, accessors and diagnostics. SC-006: reference-level projection contract
and workflow qualification. SC-007: optional compatible source-stability import.
Cross-dataset, missing/extra-program and held-out validation remain necessary.

The immediate design dependency is the benchmark review's comparison contract:
explicit units/support, valid sign orientation and a clear distinction between
loading representability, observed usage and biological program recovery.
The independent adversarial review found a prototype sign/pole defect; the
dedicated reference API should not inherit that behavior.

## Stale and current evidence

The 2026-10-16 GBCD/STARCAT backlog memo predates the shipped GBCD adapter.
Its statement that the full GBCD work is unstarted is stale for the delivered
adapter. Some older GB-005 checkboxes are stale too; those activity summaries
and display ordering already shipped. This does not imply StarCAT also shipped.

Checked current NAMESPACE, R/programs.R, R/projection.R, R/gbcd.R, HANDOFF.md,
flashier-starcat/03_EXECUTION_TRACKER.md, development/major-improvements-progress.md
and validation/gbcd-integration-20261023/REPORT.md.
No new StarCAT implementation, release, publication or scheduler work occurred.

Attribution: StarCAT-related source and method credits remain governed by the
audited design bundle and inst/provenance.csv; no unresolved individual authors
were added here. The shipped GBCD adapter credits Yusha Liu, Peter Carbonetto,
Jason Willwerscheid and Matthew Stephens, with adapter work by Mikhael Manurung.
