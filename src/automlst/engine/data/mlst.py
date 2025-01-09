from dataclasses import dataclass
from typing import Mapping, Sequence

@dataclass(frozen=True)
class Allele:
    allele_loci: str
    allele_variant: str

@dataclass(frozen=True)
class MLSTProfile:
    alleles: Mapping[str, Sequence[Allele]]
    sequence_type: str
    clonal_complex: str