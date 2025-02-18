from typing import AsyncIterable, Iterable

import pytest
from autobigs.engine.structures.alignment import AlignmentStats
from autobigs.engine.writing import alleles_to_text_map, write_mlst_profiles_as_csv
from autobigs.engine.structures.mlst import Allele, MLSTProfile, NamedMLSTProfile
import tempfile
from csv import reader
from os import path


@pytest.fixture
def dummy_alphabet_mlst_profile():
    return NamedMLSTProfile("name", MLSTProfile((
        Allele("A", "1", None),
        Allele("D", "1", None),
        Allele("B", "1", None),
        Allele("C", "1", None),
        Allele("C", "2", AlignmentStats(90, 10, 0, 90))
    ), "mysterious", "very mysterious"))

async def iterable_to_asynciterable(iterable: Iterable):
    for iterated in iterable:
        yield iterated

async def test_column_order_is_same_as_expected_file(dummy_alphabet_mlst_profile: MLSTProfile):
    dummy_profiles = [dummy_alphabet_mlst_profile]
    with tempfile.TemporaryDirectory() as temp_dir:
        output_path = path.join(temp_dir, "out.csv")
        await write_mlst_profiles_as_csv(iterable_to_asynciterable(dummy_profiles), output_path)
        with open(output_path) as csv_handle:
            csv_reader = reader(csv_handle)
            lines = list(csv_reader)
            target_columns = lines[4:]
            assert target_columns == sorted(target_columns)

async def test_alleles_to_text_map_mapping_is_correct(dummy_alphabet_mlst_profile: NamedMLSTProfile):
    mapping = alleles_to_text_map(dummy_alphabet_mlst_profile.mlst_profile.alleles) # type: ignore
    expected_mapping = {
        "A": "1",
        "B": "1",
        "C": ("1", "2*"),
        "D": "1"
    }
    for allele_name, allele_ids in mapping.items():
        assert allele_name in expected_mapping
        assert allele_ids == expected_mapping[allele_name]
