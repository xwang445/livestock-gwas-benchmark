# Submission figure scripts, 12 November 2026

These scripts reproduce the current artwork from committed non-identifying aggregates. Run from the repository root after installing `environment/requirements-submission.txt`:

```bash
python scripts/submission_2026-11-12/plot_figures_1_2_7.py --output /tmp/submission-figures
python scripts/submission_2026-11-12/plot_figures_3_5_6.py --output /tmp/submission-figures
python scripts/submission_2026-11-12/plot_figure8.py --output /tmp/submission-figures
python scripts/submission_2026-11-12/make_supplementary_table_s1.py --output /tmp/submission-figures/supplement/tables
```

The first two plotting commands accept `--figure` to select one main figure. Figure 5 also generates Figures S11–S20. The default output is `figures/submission_2026-11-12/`; use a new output directory to preserve the retained November 12 artwork. `--source` overrides the committed aggregate directory. Table S1 accepts `--tables` for the public table source directory.

Figure 8 uses final residual-variance-estimating fits for both pipelines and author-supplied **recovered-effect PVE**, rather than GWAS-coefficient PVE or certified generating effects. AP is read from exact full-ranking calculations; display vertices are used only for plotting. Figure 7 redraws retained theory bins and observations; it does not regenerate the unavailable upstream theoretical variant cache.

## Restricted-input Manhattan redraws

Figure 4 and S1–S10 require authorized variant-level association data, excluded from Git. Use the original benchmark layout and completed GCTA outputs:

```bash
python scripts/submission_2026-11-12/plot_manhattan.py   --benchmark "$GWAS_BENCHMARK_ROOT" --gcta-directory "$GCTA_OUTPUT_DIRECTORY"   --output /tmp/submission-manhattan
```

Alternatively pass `--source-cache "$AUTHORIZED_MANHATTAN_CACHE"` for the 100 local `phenoXX_METHOD.tsv.gz` tables retained with the November 12 source package. `--main-only`, `--supplement-only`, or `--phenotype 1` select outputs. Every panel is checked against the committed 100-row count table. REGENIE is tested directly on its native LOG10P scale; other methods use inclusive P ≤ 5 × 10⁻⁸. All panels use the shared 0–110 ordinate. The script writes aggregate audits with relative source filenames and hashes; it does not export variant rows.

## Export checks

`verify_submission_figures.py --root /tmp/submission-figures` checks all eight main and twenty supplementary figures and requires Poppler (`pdfinfo`, `pdffonts`, `pdftotext`, `pdftoppm`) plus Pillow. `inspect_saved_pdf_lines.py` inspects retained public PDFs for physical stroke widths. These checks require the complete figure set; aggregate-only redraws omit Figure 4 and S1–S10 until restricted inputs are supplied. Artwork exporters use Matplotlib; table DOCX generation uses the Python standard library and retains a single editable table. The original workstation scripts that edited a private manuscript are intentionally separate from this public reproduction workflow.
