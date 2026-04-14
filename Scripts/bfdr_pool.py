"""Global BFDR recomputation for merged per-group SAINT outputs.

SAINT's own BFDR is computed globally across all (Bait, Prey) pairs in one run
(verified in SAINTexpress-custom/SAINT-MRF-int/main.hpp:491-499, which iterates
Fastmat.mat flat via average_score.size()).  When we fan SAINT out per group,
each per-group list.txt carries a BFDR calibrated against only that group's
pairs.  Re-pooling restores the global posterior-FDR calibration across the
merged network.
"""

import numpy as np
import pandas as pd


def recompute_bfdr(df):
    """Recompute BFDR globally across all rows using SAINT's formula.

    For each row with AvgP = p:
        count = number of rows with AvgP > p (strict)
        BFDR  = 1 - (sum of AvgP among rows with AvgP > p) / count
        BFDR  = 0 if count == 0 (highest-scoring rows get BFDR = 0)

    Rows with identical AvgP get identical BFDR — tied rows do not contribute
    to each other's count or sum.

    Replaces the BFDR column in place.  Preserves input row order.  Leaves a
    BFDR column as float64.  Requires an 'AvgP' column; raises ValueError
    otherwise.
    """
    if "AvgP" not in df.columns:
        raise ValueError("recompute_bfdr requires an 'AvgP' column")

    out = df.copy()
    s = out["AvgP"].to_numpy(dtype=float)
    n = len(s)

    if n == 0:
        out["BFDR"] = np.array([], dtype=float)
        return out

    # Descending order by AvgP (stable so ties keep relative input order).
    order = np.argsort(-s, kind="stable")
    sorted_s = s[order]

    # For each element, count how many rows have strictly larger AvgP.
    # -sorted_s is ascending; searchsorted(..., side='left') finds the first
    # position whose -value >= -v (i.e. value <= v), which equals the number
    # of rows with value > v.  Tied rows share the same result.
    strict_greater_count = np.searchsorted(-sorted_s, -sorted_s, side="left")

    # Prefix sum of sorted values lets us read the sum of the top-k in O(1).
    prefix_sum = np.concatenate(([0.0], np.cumsum(sorted_s)))
    strict_greater_sum = prefix_sum[strict_greater_count]

    with np.errstate(divide="ignore", invalid="ignore"):
        bfdr_sorted = np.where(
            strict_greater_count > 0,
            1.0 - strict_greater_sum / strict_greater_count,
            0.0,
        )
    bfdr_sorted = np.clip(bfdr_sorted, 0.0, 1.0)

    bfdr = np.empty(n, dtype=float)
    bfdr[order] = bfdr_sorted
    out["BFDR"] = bfdr
    return out
