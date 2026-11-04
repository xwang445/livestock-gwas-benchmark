# SYNTHETIC EXAMPLE FOR CODE VALIDATION ONLY

This generator creates 32 artificial animals, 36 SNPs on two chromosomes, two phenotypes, and nine method-format fixtures. It includes effect/null windows and SNPs at exact boundaries. The format fixtures encode artificial statistics; **they are not executions of nine GWAS algorithms** and do not reproduce manuscript results.

From the repository root:

```bash
python example/synthetic/generate.py --out /tmp/synthetic-gwas-input
python -m unittest discover -s tests -v
```

Use an empty/new output directory. The unit tests generate their own fixtures and exercise all nine harmonizers, loci, spillover, lambda, overlap, power/FPR, aggregate schemas and resumable caches. Fixtures are generated deterministically with NumPy seed 20261104. Synthetic animal IDs are explicitly prefixed SYNTHETIC_ANIMAL. Real animal IDs must never be placed here.
