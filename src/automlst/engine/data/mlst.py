from dataclasses import dataclass
from typing import Mapping, Sequence, Union

@dataclass(frozen=True)
class PartialAllelicMatchProfile:
    percent_identity: float
    mismatches: int
    bitscore: float
    gaps: int

@dataclass(frozen=True)
class Allele:
    allele_loci: str
    allele_variant: str
    partial_match_profile: Union[None, PartialAllelicMatchProfile]

@dataclass(frozen=True)
class MLSTProfile:
    alleles: Mapping[str, Sequence[Allele]]
    sequence_type: str
    clonal_complex: str
