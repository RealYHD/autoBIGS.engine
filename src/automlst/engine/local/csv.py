import csv
from io import TextIOWrapper
from os import PathLike
from typing import AsyncIterable, Iterable, Mapping, Sequence, Union

from automlst.engine.data.mlst import Allele, MLSTProfile


def dict_loci_alleles_variants_from_loci(alleles_map: Mapping[str, Sequence[Allele]]):
    result_dict: dict[str, list[str]] = {}
    for loci, alleles in alleles_map.items():
        result_dict[loci] = list()
        for allele in alleles:
            result_dict[loci].append(allele.allele_variant)
    return result_dict


async def write_mlst_profiles_as_csv(mlst_profiles_iterable: AsyncIterable[tuple[str, Union[MLSTProfile, None]]], handle: Union[str, bytes, PathLike[str], PathLike[bytes]]) -> Sequence[str]:
    failed = list()
    with open(handle, "w", newline='') as filehandle:
        header = None
        writer: Union[csv.DictWriter, None] = None
        async for name, mlst_profile in mlst_profiles_iterable:
            if mlst_profile is None:
                failed.append(name)
                continue
            if writer is None:
                header = ["id", "st", "clonal-complex", *mlst_profile.alleles.keys()]
                writer = csv.DictWriter(filehandle, fieldnames=header)
                writer.writeheader()
            row_dictionary = {
                "st": mlst_profile.sequence_type,
                "clonal-complex": mlst_profile.clonal_complex,
                "id": name,
                **dict_loci_alleles_variants_from_loci(mlst_profile.alleles)
            }
            writer.writerow(rowdict=row_dictionary)
    return failed