import csv
from io import TextIOWrapper
from os import PathLike
from typing import AsyncIterable, Iterable, Mapping, Sequence, Union

from automlst.engine.data.MLST import Allele, MLSTProfile


def loci_alleles_variants_from_loci(alleles_map: Mapping[str, Sequence[Allele]]):
    result_dict: dict[str, list[str]] = {}
    for loci, alleles in alleles_map.items():
        result_dict[loci] = list()
        for allele in alleles:
            result_dict[loci].append(allele.allele_variant)
    return result_dict


async def write_mlst_profiles_as_csv(mlst_profiles_iterable: Iterable[MLSTProfile], handle: Union[str, bytes, PathLike[str], PathLike[bytes]]):
    mlst_profiles = list(mlst_profiles_iterable)
    header = ["st", "clonal-complex", *mlst_profiles[0].alleles.keys()]
    with open(handle, "w", newline='') as filehandle:
        writer = csv.DictWriter(filehandle, fieldnames=header)
        writer.writeheader()
        for mlst_profile in mlst_profiles:
            row_dictionary = {
                "st": mlst_profile.sequence_type,
                "clonal-complex": mlst_profile.clonal_complex,
                **loci_alleles_variants_from_loci(mlst_profile.alleles)
            }

            writer.writerow(rowdict=row_dictionary)
