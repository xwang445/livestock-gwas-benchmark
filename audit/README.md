# Reconstruction audit

Audit date: 2026-11-04. The selected manuscript is the user-confirmed October 22 DOCX, snapshot SHA256 `691b4cad7f5a771392009a949becd0c8dded2b99a6481bb47f8069bc2dc34e09` (modified November 4). Its full text, equations and representative rendered figure pages were inspected. The July bioRxiv version does not include the current fine-mapping section.

Local starting commit: `ec6cb3297cbcb1e1756da1a76a1f7be8ad626088`, clean working tree. Public baseline fetched read-only: `14ff11e8748f2fbca22926e555edf4763942c9f2` (October 20 fine-mapping update). Work is on a separate local review branch. No push.

Full metadata inventories are local audit records at `audit/inventory_local.tsv`, `audit/inventory_hpcc.tsv`, and `audit/inventory_manuscript.tsv`. They are deliberately Git-ignored: they contain restricted-source paths and over 500,000 HPCC entries. `inventory_summary.tsv` records counts and SHA256 fingerprints. The public Git inventory and code provenance table are reviewable here. Source JSON, unmodified code snapshots, manuscript snapshot and renders are retained in the private audit working directory recorded in local PROJECT_STATE.md; they are not publication data.

Classification is evidence-based: a candidate is not authoritative merely because it is newer. UNKNOWN remains explicit for unverified candidates. Prior audit reports are evidence to investigate, not a substitute for independent validation. Raw genotypes, animal IDs, phenotypes, GRMs, LD matrices, per-variant GWAS and fine-mapping calls are excluded from tracked outputs.
