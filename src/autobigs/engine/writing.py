from collections import defaultdict
import csv
from os import PathLike
from typing import AsyncIterable, Collection, Mapping, Sequence, Union

from autobigs.engine.structures.mlst import Allele, MLSTProfile, NamedMLSTProfile


def alleles_to_text_map(alleles: Collection[Allele]) -> Mapping[str, Union[Sequence[str], str]]:
    result = defaultdict(list)
    for allele in alleles:
        result[allele.allele_locus].append(allele.allele_variant + ("*" if allele.partial_match_profile is not None else ""))
    for locus in result.keys():
        if len(result[locus]) == 1:
            result[locus] = result[locus][0] # Take the only one
        else:
            result[locus] = tuple(result[locus]) # type: ignore
    return dict(result)

async def write_mlst_profiles_as_csv(mlst_profiles_iterable: AsyncIterable[NamedMLSTProfile], handle: Union[str, bytes, PathLike[str], PathLike[bytes]]) -> Sequence[str]:
    failed = list()
    with open(handle, "w", newline='') as filehandle:
        header = None
        writer: Union[csv.DictWriter, None] = None
        async for named_mlst_profile in mlst_profiles_iterable:
            name = named_mlst_profile.name
            mlst_profile = named_mlst_profile.mlst_profile
            if mlst_profile is None:
                failed.append(name)
                continue
            allele_mapping = alleles_to_text_map(mlst_profile.alleles)
            if writer is None:
                header = ["id", "st", "clonal-complex", *sorted(allele_mapping.keys())]
                writer = csv.DictWriter(filehandle, fieldnames=header)
                writer.writeheader()
            row_dictionary = {
                "st": mlst_profile.sequence_type,
                "clonal-complex": mlst_profile.clonal_complex,
                "id": name,
                **allele_mapping
            }
            writer.writerow(rowdict=row_dictionary)
    return failed