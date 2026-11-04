# Decisions

- 2026-11-04: use user-confirmed October22 manuscript, not July bioRxiv draft.
- Preserve original local main and use a review branch based on October20 public HEAD.
- Full path inventories remain local and ignored; publish only code provenance and aggregate counts.
- No fresh GWAS or fine-mapping fitting: recover and independently recompute downstream summaries first.
- Preserve strict Figure3 versus inclusive spillover thresholds as separate documented historical definitions.
- Use evidence-supported GAPIT Kansas/NYC branches; do not infer from names or modify p-values.
- Do not report the recovered PVE estimator as an original sampled effect or force it to match the manuscript.
- Do not call downstream Equation13 checks an end-to-end validation of missing Equation12 code.

- Independent reconstruction matched strict Figure 3, historical lambda, power/FPR/spatial, Equation 13 downstream and fine-mapping metrics; retain original rules and disclose their limitations.
- Current GAPIT source version does not prove executed version. Pin retained bytes for future reviewable runs, and keep historical version unresolved.
- Two numerical manuscript discrepancies remain FAIL: REGENIE maximum count and full-range Manhattan maximum. No result editing to force agreement.
- Rebuilt outputs after adding source-tree write guards; original source folders remain untouched.
- Added local-review figures and aggregates only; any eventual push requires separate human review/authorization and data-rights confirmation.

- 2026-11-04 publication update: after the completion report, the user explicitly authorized pushing the reviewed package to GitHub. Pushed the five reviewed commits through d49a978 to origin/reconstruct/20261104. No raw restricted data were included. Scientific unresolved items remain unchanged.

- 2026-11-04 default-branch correction: the user found only templates because the prior authorized push updated the review branch, not main. Publish the reviewed package to main, preserve all five templates under example/schemas, and make derived_data contain executed results only. Original data and numeric results remain unchanged.

- 2026-11-04 configuration update: generated the result snapshot from executed tables, verified all 90 Mac/HPCC GWAS paths, and recorded machine-specific input availability. Absolute roots/SSH host are kept in ignored local_paths.yaml; public scripts consume relative paths. Unknown historical settings and original simulation identity remain unresolved.
