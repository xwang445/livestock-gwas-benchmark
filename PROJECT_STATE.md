# Project state

Status: initial audit documented; independent aggregate validation and implementation next. Not complete. No push.

Baseline: local clean ec6cb32; public14ff11e. Review branch reconstruct/20261104. Original Mac and HPCC sources remain unchanged. Private audit workspace: /private/tmp/livestock-reconstruction-20261104.

Completed: environment/SSH checks; recursive inventories (590 local,502804HPCC,69public files); upstream source/log selection; user-confirmed October22 manuscript snapshot and full text/equation inspection; representative figure rendering; candidate provenance and discrepancy ledger.

Active jobs: none. No expensive GWAS or fine-mapping fits submitted.

Failures: git switch unsupported by Git2.22; recovered using checkout. fitz absent from document runtime; used bundled pypdf and official renderer successfully. Initial sandbox network denial resolved by approved read-only network commands.

Open: upstream Equation12 generator, exact maximum-PVE claim provenance, array metadata and production imputation version; current source order discrepancy; stage-specific significance/interval definitions.

Next: implement scripts/benchmark_core.py, run all90existing local GWAS through independent downstream calculations into a new output directory, then compare against saved manuscript tables. Resume using source hashes/cache, never restart completed jobs blindly.
