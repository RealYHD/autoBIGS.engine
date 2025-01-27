from dataclasses import dataclass
from typing import Mapping, Sequence, Union

@dataclass(frozen=True)
class PartialAllelicMatchProfile:
    percent_identity: float
    mismatches: int
    gaps: int

@dataclass(frozen=True)
class Allele:
    allele_locus: str
    allele_variant: str
    partial_match_profile: Union[None, PartialAllelicMatchProfile]

@dataclass(frozen=True)
class MLSTProfile:
    alleles: Mapping[str, Allele]
    sequence_type: str
    clonal_complex: str

@dataclass(frozen=True)
class NamedMLSTProfile:
    name: str
    mlst_profile: Union[None, MLSTProfile]