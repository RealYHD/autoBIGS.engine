from typing import AsyncIterable, Iterable
from autobigs.engine.data.local.csv import dict_loci_alleles_variants_from_loci, write_mlst_profiles_as_csv
from autobigs.engine.data.structures.mlst import Allele, MLSTProfile
import tempfile
from csv import reader
from os import path

async def iterable_to_asynciterable(iterable: Iterable):
    for iterated in iterable:
        yield iterated

def test_dict_loci_alleles_variants_from_loci_single_loci_not_list():
    alleles_map = {
        "adk": [Allele("adk", "1", None)]
    }
    results = dict_loci_alleles_variants_from_loci(alleles_map)
    for loci, variant in results.items():
        assert isinstance(variant, str)
        assert variant == "1"

def test_dict_loci_alleles_variants_from_loci_multi_loci_is_list():
    alleles_map = {
        "adk": [Allele("adk", "1", None), Allele("adk", "2", None)]
    }
    results = dict_loci_alleles_variants_from_loci(alleles_map)
    for loci, variant in results.items():
        assert isinstance(variant, list)
        assert len(variant) == 2

async def test_column_order_is_same_as_expected_file():
    dummy_profiles = [("test_1", MLSTProfile({
        "A": Allele("A", "1", None),
        "D": Allele("D", "1", None),
        "B": Allele("B", "1", None),
        "C": Allele("C", "1", None)
    }, "mysterious", "very mysterious"))]
    with tempfile.TemporaryDirectory() as temp_dir:
        output_path = path.join(temp_dir, "out.csv")
        await write_mlst_profiles_as_csv(iterable_to_asynciterable(dummy_profiles), output_path)
        with open(output_path) as csv_handle:
            csv_reader = reader(csv_handle)
            lines = list(csv_reader)
            target_columns = lines[4:]
            assert target_columns == sorted(target_columns)            