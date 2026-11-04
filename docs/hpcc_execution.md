# HPCC execution

The supplied analysis host was reachable. R, Python, PLINK and PLINK2 were present; `sbatch`, `srun` and `qsub` were not on its PATH. No scheduler, account, partition or authorized resource allocation can be inferred from that. No expensive fitting jobs were submitted.

The audit used read-only commands on the historical tree. Independent LD verification streamed the existing compressed pairwise LD file with one BLAS/OMP thread; fine-mapping verification read saved fit summaries/PIPs and MAF annotations. Only aggregate JSON was returned. Local downstream summaries used saved local GWAS outputs. All these verification jobs completed.

For future fitting, create a separate authorized workspace. Follow scripts/12_finemapping/README.md for installer and dependencies. Genotype preparation precedes simulation; simulation and sample-identity checks precede GWAS; harmonization and summary statistics follow GWAS; LOCO adjusted LD precedes SuSiE-adj; ordinary LD precedes RSS. Truth labels enter evaluation only.

Historical GWAS wrappers used 16 threads (MPH simulation 20). LD logging shows up to 79 threads and automatic memory reservation; the recovered LD wrapper instead requires explicit approved THREADS and MEMORY_MB to prevent uncontrolled reservation. Fine-mapping includes dense 20,000 × 20,000 matrices: 3.2 GB per float64 matrix before extra work arrays. Existing fit runners use substantial concurrency (LOCO 40 R workers; RSS 20); set concurrency to the available allocation before execution. No memory/runtime claim has been benchmarked in this reconstruction.

The SLURM template below is **prepared, not submitted**, and is only applicable after the actual scheduler is confirmed. The resource values and account/partition must be supplied by the authorized operator; they are not historical scientific parameters.

```bash
sbatch --cpus-per-task=APPROVED_CPUS --mem=APPROVED_MEMORY --time=APPROVED_TIME --export=ALL scripts/hpcc/run_stage.sbatch /path/to/repository/scripts/04_gwas/run_boltlmm.sh
```

Exact stage commands and environment inputs are documented in workflow.md and each runner. Capture software versions, scheduler job IDs, source hashes, elapsed time, exit status and output locations in audit/run_registry.tsv for any future execution. Never install a new run over historical fits.
