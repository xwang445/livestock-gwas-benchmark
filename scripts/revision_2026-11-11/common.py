"""Shared validation and path-free provenance for recovered-effect analyses."""
import hashlib
from pathlib import Path

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(role, path):
    """Record a source's role and content, without its private pathname."""
    return dict(role=role, bytes=Path(path).stat().st_size, sha256=sha256(path))


def causal_ids(path):
    import pandas as pd
    truth = pd.read_csv(path)
    require('SNP' in truth, 'Causal CSV requires a SNP column.')
    require(truth.SNP.notna().all() and truth.SNP.is_unique, 'Causal SNP IDs must be unique.')
    return set(truth.SNP)


def effect_long(effect, causal, n_phenotypes=10):
    """Retain supplied PVE values; do not infer original generating effects."""
    import pandas as pd
    columns = ['pve{}'.format(k) for k in range(1, n_phenotypes + 1)]
    require(effect.snp.is_unique and effect.snp.notna().all(), 'Repeated or missing recovered SNP IDs.')
    require(set(effect.snp) == causal, 'Recovered SNPs must equal the fixed causal truth.')
    values = effect[columns].to_numpy(float)
    require(np.isfinite(values).all() and (values >= 0).all(), 'Recovered PVE must be finite and nonnegative.')
    return pd.concat([
        pd.DataFrame(dict(SNP=effect.snp, phenotype=k,
                          recovered_pve_fraction=effect['pve{}'.format(k)]))
        for k in range(1, n_phenotypes + 1)
    ], ignore_index=True)
