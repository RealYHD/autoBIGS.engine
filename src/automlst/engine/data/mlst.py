from dataclasses import dataclass
from typing import Mapping, Sequence

@dataclass
class Allele:
    allele_loci: str
    allele_variant: str

@dataclass
class MLSTProfile:
    alleles: Mapping[str, Sequence[Allele]]
    sequence_type: int
    clonal_complex: str