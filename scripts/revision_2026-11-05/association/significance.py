"""One association decision rule for the current revision.

Compare in the source scale. In particular, never exponentiate native -log10(P)
values before applying the boundary, and never widen it with an isclose tolerance.
Strict comparisons exist only for explicitly labeled historical audit columns.
"""
import numpy as np

ALPHA = 5e-8
COMPARATOR = 'le'
RULE = 'P <= 5e-8'


def significant(values, is_log10=False, alpha=ALPHA):
    """Finite valid raw P <= alpha, or native -log10(P) >= -log10(alpha)."""
    alpha=float(alpha)
    if not 0 < alpha <= 1:
        raise ValueError('Association alpha must be in (0, 1]')
    v=np.asarray(values,dtype=float)
    valid=np.isfinite(v)&(v>=0)
    if is_log10:
        return valid&(v>=-np.log10(alpha))
    return valid&(v<=1)&(v<=alpha)


def strict_for_historical_audit(values, is_log10=False, alpha=ALPHA):
    """Legacy comparison only; never supplies current primary result columns."""
    v=np.asarray(values,dtype=float)
    valid=np.isfinite(v)&(v>=0)
    if is_log10:
        return valid&(v>-np.log10(alpha))
    return valid&(v<=1)&(v<alpha)
