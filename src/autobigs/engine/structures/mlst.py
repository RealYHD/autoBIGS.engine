from collections import defaultdict
from dataclasses import dataclass
from typing import Collection, Iterable, Mapping, Sequence, Union

from autobigs.engine.structures.alignment import AlignmentStats

@dataclass(frozen=True)
class Allele:
    allele_locus: str
    allele_variant: str
    partial_match_profile: Union[None, AlignmentStats]

@dataclass(frozen=True)
class MLSTProfile:
    alleles: Collection[Allele]
    sequence_type: str
    clonal_complex: str

@dataclass(frozen=True)
class NamedMLSTProfile:
    name: str
    mlst_profile: Union[None, MLSTProfile]


def alleles_to_mapping(alleles: Iterable[Allele]):
    result = defaultdict(list)
    for allele in alleles:
        result[allele.allele_locus].append(allele.allele_variant)
    result = dict(result)
    for locus, variant in result.items():
        if len(variant) == 1:
            result[locus] = variant[0]
    return result