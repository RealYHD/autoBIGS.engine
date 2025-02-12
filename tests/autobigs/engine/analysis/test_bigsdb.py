from os import path
import random
import re
from typing import Callable, Collection, Sequence, Union
from Bio import SeqIO
import pytest
from autobigs.engine.analysis import bigsdb
from autobigs.engine.structures import mlst
from autobigs.engine.structures.genomics import NamedString
from autobigs.engine.structures.mlst import Allele, MLSTProfile
from autobigs.engine.exceptions.database import NoBIGSdbExactMatchesException, NoBIGSdbMatchesException
from autobigs.engine.analysis.bigsdb import BIGSdbIndex, BIGSdbMLSTProfiler, RemoteBIGSdbMLSTProfiler

async def generate_async_iterable(normal_iterable):
    for dummy_sequence in normal_iterable:
        yield dummy_sequence

def gene_scrambler(gene: str, mutation_site_count: Union[int, float], alphabet: Sequence[str] = ["A", "T", "C", "G"]):
    rand = random.Random(gene)
    if isinstance(mutation_site_count, float):
        mutation_site_count = int(mutation_site_count * len(gene))
    random_locations = rand.choices(range(len(gene)), k=mutation_site_count)
    scrambled = list(gene)
    for random_location in random_locations:
        scrambled[random_location] = rand.choice(alphabet)
    return "".join(scrambled)

def get_first_sequence_from_fasta(resource: str):
    return str(SeqIO.read(path.join("tests/resources/", resource), "fasta").seq)

def get_multiple_sequences_from_fasta(resource: str):
    return tuple(SeqIO.parse(path.join("tests/resources/", resource), "fasta"))

bpertussis_tohamaI_profile = MLSTProfile((
        Allele("adk", "1", None),
        Allele("fumC", "1", None),
        Allele("glyA", "1", None),
        Allele("tyrB", "1", None),
        Allele("icd", "1", None),
        Allele("pepA", "1", None),
        Allele("pgm", "1", None)), "1", "ST-2 complex")

bpertussis_tohamaI_bad_profile = MLSTProfile((
        Allele("adk", "1", None),
        Allele("fumC", "2", None),
        Allele("glyA", "36", None),
        Allele("tyrB", "4", None),
        Allele("icd", "4", None),
        Allele("pepA", "1", None),
        Allele("pgm", "5", None),
    ), "unknown", "unknown")

hinfluenzae_2014_102_profile = MLSTProfile((
        Allele("adk", "28", None),
        Allele("atpG", "33", None),
        Allele("frdB", "7", None),
        Allele("fucK", "18", None),
        Allele("mdh", "11", None),
        Allele("pgi", "125", None),
        Allele("recA", "89", None)
    ), "478", "unknown")

hinfluenzae_2014_102_bad_profile = MLSTProfile((
        Allele("adk", "3", None),
        Allele("atpG", "121", None),
        Allele("frdB", "6", None),
        Allele("fucK", "5", None),
        Allele("mdh", "12", None),
        Allele("pgi", "4", None),
        Allele("recA", "5", None)
    ), "unknown", "unknown")


@pytest.mark.parametrize("local_db,database_api,database_name,schema_id,seq_path,feature_seqs_path,expected_profile,bad_profile", [
    (False, "https://bigsdb.pasteur.fr/api", "pubmlst_bordetella_seqdef", 3, "tohama_I_bpertussis.fasta", "tohama_I_bpertussis_features.fasta", bpertussis_tohamaI_profile, bpertussis_tohamaI_bad_profile),
    (False, "https://rest.pubmlst.org", "pubmlst_hinfluenzae_seqdef", 1, "2014-102_hinfluenza.fasta", "2014-102_hinfluenza_features.fasta", hinfluenzae_2014_102_profile, hinfluenzae_2014_102_bad_profile),
])
class TestBIGSdbMLSTProfiler:
    async def test_profiling_results_in_exact_matches_when_exact(self, local_db, database_api, database_name, schema_id, seq_path: str, feature_seqs_path: str, expected_profile: MLSTProfile, bad_profile: MLSTProfile):
        sequence = get_first_sequence_from_fasta(seq_path)
        async with bigsdb.get_BIGSdb_MLST_profiler(local_db, database_api, database_name, schema_id) as dummy_profiler:
            expected_alleles = mlst.alleles_to_mapping(expected_profile.alleles)
            targets_left = set(mlst.alleles_to_mapping(expected_profile.alleles).keys())
            async for exact_match in dummy_profiler.determine_mlst_allele_variants(query_sequence_strings=[sequence]):
                assert isinstance(exact_match, Allele)
                assert exact_match.allele_locus in expected_alleles
                assert exact_match.allele_variant == expected_alleles[exact_match.allele_locus]
                targets_left.remove(exact_match.allele_locus)

            assert len(targets_left) == 0

    async def test_sequence_profiling_non_exact_returns_non_exact(self, local_db, database_api, database_name, schema_id, seq_path: str, feature_seqs_path: str, expected_profile: MLSTProfile, bad_profile: MLSTProfile):
        target_sequences = get_multiple_sequences_from_fasta(feature_seqs_path)
        mlst_targets = {x.lower() for x in mlst.alleles_to_mapping(expected_profile.alleles).keys()}
        async with bigsdb.get_BIGSdb_MLST_profiler(local_db, database_api, database_name, schema_id) as profiler:
            for target_sequence in target_sequences:
                match = re.fullmatch(r".*\[gene=([\w\d]+)\].*", target_sequence.description)
                if match is None:
                    continue
                gene = match.group(1).lower()
                if gene not in mlst_targets:
                    continue
                scrambled = gene_scrambler(str(target_sequence.seq), 0.125)
                async for partial_match in profiler.determine_mlst_allele_variants([scrambled]):
                    assert partial_match.partial_match_profile is not None
                    mlst_targets.remove(gene)

            assert len(mlst_targets) == 0

    async def test_profiling_results_in_correct_mlst_st(self, local_db, database_api, database_name, schema_id, seq_path: str, feature_seqs_path: str, expected_profile: MLSTProfile, bad_profile: MLSTProfile):
        async with bigsdb.get_BIGSdb_MLST_profiler(local_db, database_api, database_name, schema_id) as dummy_profiler:
            mlst_st_data = await dummy_profiler.determine_mlst_st(expected_profile.alleles)
            assert mlst_st_data is not None
            assert isinstance(mlst_st_data, MLSTProfile)
            assert mlst_st_data.clonal_complex == expected_profile.clonal_complex
            assert mlst_st_data.sequence_type == expected_profile.sequence_type

    async def test_profiling_non_exact_results_in_list_of_mlsts(self, local_db, database_api, database_name, schema_id, seq_path: str, feature_seqs_path: str, expected_profile: MLSTProfile, bad_profile: MLSTProfile):
        dummy_alleles = bad_profile.alleles
        async with bigsdb.get_BIGSdb_MLST_profiler(local_db, database_api, database_name, schema_id) as dummy_profiler:
            mlst_profile = await dummy_profiler.determine_mlst_st(dummy_alleles)
            assert mlst_profile.clonal_complex == "unknown"
            assert mlst_profile.sequence_type == "unknown"


    async def test_bigsdb_profile_multiple_strings_same_string_twice(self, local_db, database_api, database_name, schema_id, seq_path: str, feature_seqs_path: str, expected_profile: MLSTProfile, bad_profile: MLSTProfile):
        sequence = get_first_sequence_from_fasta(seq_path)
        dummy_sequences = [[NamedString("seq1", sequence)], [NamedString("seq2", sequence)]]
        async with bigsdb.get_BIGSdb_MLST_profiler(local_db, database_api, database_name, schema_id) as dummy_profiler:
            async for named_profile in dummy_profiler.profile_multiple_strings(generate_async_iterable(dummy_sequences)):
                name, profile = named_profile.name, named_profile.mlst_profile
                assert profile is not None
                assert isinstance(profile, MLSTProfile)
                assert profile.clonal_complex == expected_profile.clonal_complex
                assert profile.sequence_type == expected_profile.sequence_type

    async def test_bigsdb_profile_multiple_strings_exactmatch_fail_second_no_stop(self, local_db, database_api, database_name, schema_id, seq_path: str, feature_seqs_path: str, expected_profile: MLSTProfile, bad_profile: MLSTProfile):
        valid_seq = get_first_sequence_from_fasta(seq_path)
        dummy_sequences = [[NamedString("seq1", valid_seq)], [NamedString("should_fail", gene_scrambler(valid_seq, 0.3))], [NamedString("seq3", valid_seq)]]
        async with bigsdb.get_BIGSdb_MLST_profiler(local_db, database_api, database_name, schema_id) as dummy_profiler:
            async for name_profile in dummy_profiler.profile_multiple_strings(generate_async_iterable(dummy_sequences), True):
                name, profile = name_profile.name, name_profile.mlst_profile

                assert profile is not None
                if name == "should_fail":
                    assert profile.clonal_complex == "unknown"
                    assert profile.sequence_type == "unknown"
                    assert len(profile.alleles) > 0
                else:
                    assert isinstance(profile, MLSTProfile)
                    assert profile.clonal_complex == expected_profile.clonal_complex
                    assert profile.sequence_type == expected_profile.sequence_type

    async def test_bigsdb_profile_multiple_strings_nonexact_second_no_stop(self, local_db, database_api, database_name, schema_id, seq_path: str, feature_seqs_path: str, expected_profile: MLSTProfile, bad_profile: MLSTProfile):
        valid_seq = get_first_sequence_from_fasta(seq_path)
        dummy_sequences = [[NamedString("seq1", valid_seq)], [NamedString("should_fail", gene_scrambler(valid_seq, 0.3))], [NamedString("seq3", valid_seq)]]

        async with bigsdb.get_BIGSdb_MLST_profiler(local_db, database_api, database_name, schema_id) as dummy_profiler:
            async for named_profile in dummy_profiler.profile_multiple_strings(generate_async_iterable(dummy_sequences), False):
                name, profile = named_profile.name, named_profile.mlst_profile
                
                assert profile is not None
                if name == "should_fail":
                    assert profile.clonal_complex == "unknown"
                    assert profile.sequence_type == "unknown"
                    assert len(profile.alleles) > 0
                else:
                    assert isinstance(profile, MLSTProfile)
                    assert profile.clonal_complex == expected_profile.clonal_complex
                    assert profile.sequence_type == expected_profile.sequence_type

class TestBIGSdbIndex:

    async def test_bigsdb_index_all_databases_is_not_empty(self):
        async with BIGSdbIndex() as bigsdb_index:
            assert len(await bigsdb_index.get_known_seqdef_dbs()) > 0

    async def test_bigsdb_index_references_pubmlst_correctly(self):
        async with BIGSdbIndex() as bigsdb_index:
            assert (await bigsdb_index.get_bigsdb_api_from_seqdefdb("pubmlst_hinfluenzae_seqdef")) == "https://rest.pubmlst.org"

    async def test_bigsdb_index_references_institutpasteur_correctly(self):
        async with BIGSdbIndex() as bigsdb_index:
            assert (await bigsdb_index.get_bigsdb_api_from_seqdefdb("pubmlst_bordetella_seqdef")) == "https://bigsdb.pasteur.fr/api"

    async def test_bigsdb_index_get_schemas_for_bordetella(self):
        async with BIGSdbIndex() as index:
            schemas = await index.get_schemas_for_seqdefdb(seqdef_db_name="pubmlst_bordetella_seqdef")
            assert len(schemas.keys()) > 0
            assert "MLST" in schemas
            assert isinstance(schemas["MLST"], int)

    async def test_bigsdb_index_get_databases_has_only_seqdef(self):
        async with BIGSdbIndex() as index:
            databases = await index.get_known_seqdef_dbs()
            assert len(databases.keys()) > 0
            for database_name in databases.keys():
                assert database_name.endswith("seqdef")
            assert databases["pubmlst_bordetella_seqdef"] == "https://bigsdb.pasteur.fr/api"

    @pytest.mark.parametrize("local", [
        (True),
        (False)
    ])
    async def test_bigsdb_index_instantiates_correct_profiler(self, local):
        sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
        async with BIGSdbIndex() as bigsdb_index:
            async with await bigsdb_index.build_profiler_from_seqdefdb(local, "pubmlst_bordetella_seqdef", 3) as profiler:
                assert isinstance(profiler, BIGSdbMLSTProfiler)
                profile = await profiler.profile_string(sequence)
                assert profile.clonal_complex == "ST-2 complex"
                assert profile.sequence_type == "1"
