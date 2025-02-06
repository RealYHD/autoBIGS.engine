from dataclasses import dataclass
from numbers import Number
from typing import Sequence

@dataclass(frozen=True)
class AlignmentStats:
    percent_identity: float
    mismatches: int
    gaps: int
    match_metric: int

@dataclass(frozen=True)
class PairwiseAlignment:
    reference: str
    query: str
    reference_indices: Sequence[Number]
    query_indices: Sequence[Number]
    alignment_stats: AlignmentStats