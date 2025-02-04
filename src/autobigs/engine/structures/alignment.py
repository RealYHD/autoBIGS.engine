from dataclasses import dataclass
from numbers import Number

@dataclass(frozen=True)
class AlignmentStats:
    percent_identity: float
    mismatches: int
    gaps: int
    score: int

@dataclass(frozen=True)
class PairwiseAlignment:
    reference: str
    query: str
    reference_indices: list[Number]
    query_indices: list[Number]
    alignment_stats: AlignmentStats