import random
import re
from typing import Collection, Sequence, Union
from Bio import SeqIO
import pytest
from autobigs.engine.data.structures.genomics import NamedString
from autobigs.engine.data.structures.mlst import Allele, MLSTProfile
from autobigs.engine.exceptions.database import NoBIGSdbExactMatchesException, NoBIGSdbMatchesException
from autobigs.engine.data.remote.databases.bigsdb import BIGSdbIndex, OnlineBIGSdbMLSTProfiler

def gene_scrambler(gene: str, mutation_site_count: Union[int, float], alphabet: Sequence[str] = ["A", "T", "C", "G"]):
    rand = random.Random(gene)
    if isinstance(mutation_site_count, float):
        mutation_site_count = int(mutation_site_count * len(gene))
    random_locations = rand.choices(range(len(gene)), k=mutation_site_count)
    scrambled = list(gene)
    for random_location in random_locations:
        scrambled[random_location] = rand.choice(alphabet)
    return "".join(scrambled)

async def test_institutpasteur_profiling_results_in_exact_matches_when_exact():
    sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    async with OnlineBIGSdbMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        targets_left = {"adk", "fumC", "glyA", "tyrB", "icd", "pepA", "pgm"}
        async for exact_match in dummy_profiler.fetch_mlst_allele_variants(sequence_strings=[sequence]):
            assert isinstance(exact_match, Allele)
            assert exact_match.allele_variant == '1' # All of Tohama I has allele id I
            targets_left.remove(exact_match.allele_locus)

        assert len(targets_left) == 0

async def test_institutpasteur_sequence_profiling_non_exact_returns_non_exact():
    sequences = list(SeqIO.parse("tests/resources/tohama_I_bpertussis_coding.fasta", "fasta"))
    mlst_targets = {"adk", "fumc", "glya", "tyrb", "icd", "pepa", "pgm"}
    async with OnlineBIGSdbMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as profiler:
        for sequence in sequences:
            match = re.fullmatch(r".*\[gene=([\w\d]+)\].*", sequence.description)
            if match is None:
                continue
            gene = match.group(1)
            if gene.lower() not in mlst_targets:
                continue
            scrambled = gene_scrambler(str(sequence.seq), 0.125)
            async for partial_match in profiler.fetch_mlst_allele_variants(scrambled):
                assert partial_match.partial_match_profile is not None
                mlst_targets.remove(gene.lower())

        assert len(mlst_targets) == 0

async def test_institutpasteur_profiling_results_in_correct_mlst_st():
    async def dummy_allele_generator():
        dummy_alleles = [
        Allele("adk", "1", None),
        Allele("fumC", "1", None),
        Allele("glyA", "1", None),
        Allele("tyrB", "1", None),
        Allele("icd", "1", None),
        Allele("pepA", "1", None),
        Allele("pgm", "1", None),
        ]
        for dummy_allele in dummy_alleles:
            yield dummy_allele
    async with OnlineBIGSdbMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        mlst_st_data = await dummy_profiler.fetch_mlst_st(dummy_allele_generator())
        assert mlst_st_data is not None
        assert isinstance(mlst_st_data, MLSTProfile)
        assert mlst_st_data.clonal_complex == "ST-2 complex"
        assert mlst_st_data.sequence_type == "1"

async def test_institutpasteur_profiling_non_exact_results_in_list_of_mlsts():
    dummy_alleles = [
    Allele("adk", "1", None),
    Allele("fumC", "2", None),
    Allele("glyA", "36", None),
    Allele("tyrB", "4", None),
    Allele("icd", "4", None),
    Allele("pepA", "1", None),
    Allele("pgm", "5", None),
    ]
    async with OnlineBIGSdbMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        mlst_profile = await dummy_profiler.fetch_mlst_st(dummy_alleles)
        assert mlst_profile.clonal_complex == "unknown"
        assert mlst_profile.sequence_type == "unknown"


async def test_institutpasteur_sequence_profiling_is_correct():
    sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    async with OnlineBIGSdbMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        profile = await dummy_profiler.profile_string(sequence)
        assert profile is not None
        assert isinstance(profile, MLSTProfile)
        assert profile.clonal_complex == "ST-2 complex"
        assert profile.sequence_type == "1"
    

async def test_pubmlst_profiling_results_in_exact_matches_when_exact():
    dummy_alleles = {
        Allele("adk", "1", None),
        Allele("atpG", "1", None),
        Allele("frdB", "1", None),
        Allele("fucK", "1", None),
        Allele("mdh", "1", None),
        Allele("pgi", "1", None),
        Allele("recA", "5", None),
    }
    sequence = str(SeqIO.read("tests/resources/FDAARGOS_1560.fasta", "fasta").seq)
    async with OnlineBIGSdbMLSTProfiler(database_api="https://rest.pubmlst.org/", database_name="pubmlst_hinfluenzae_seqdef", schema_id=1) as dummy_profiler:
        exact_matches = dummy_profiler.fetch_mlst_allele_variants(sequence_strings=sequence)
        async for exact_match in exact_matches:
            assert isinstance(exact_match, Allele)
            dummy_alleles.remove(exact_match)

        assert len(dummy_alleles) == 0

async def test_pubmlst_profiling_results_in_correct_st():
    async def generate_dummy_targets():
        dummy_alleles = [
                Allele("adk", "1", None),
                Allele("atpG", "1", None),
                Allele("frdB", "1", None),
                Allele("fucK", "1", None),
                Allele("mdh", "1", None),
                Allele("pgi", "1", None),
                Allele("recA", "5", None),
            ]
        for dummy_allele in dummy_alleles:
            yield dummy_allele
    async with OnlineBIGSdbMLSTProfiler(database_api="https://rest.pubmlst.org/", database_name="pubmlst_hinfluenzae_seqdef", schema_id=1) as dummy_profiler:
        mlst_st_data = await dummy_profiler.fetch_mlst_st(generate_dummy_targets())
        assert mlst_st_data is not None
        assert isinstance(mlst_st_data, MLSTProfile)
        assert mlst_st_data.clonal_complex == "ST-3 complex"
        assert mlst_st_data.sequence_type == "3"

async def test_pubmlst_sequence_profiling_is_correct():
    sequence = str(SeqIO.read("tests/resources/FDAARGOS_1560.fasta", "fasta").seq)
    async with OnlineBIGSdbMLSTProfiler(database_api="https://rest.pubmlst.org/", database_name="pubmlst_hinfluenzae_seqdef", schema_id=1) as dummy_profiler:
        profile = await dummy_profiler.profile_string(sequence)
        assert profile is not None
        assert isinstance(profile, MLSTProfile)
        assert profile.clonal_complex == "ST-3 complex"
        assert profile.sequence_type == "3"

async def test_bigsdb_index_all_databases_is_not_empty():
    async with BIGSdbIndex() as bigsdb_index:
        assert len(await bigsdb_index.get_known_seqdef_dbs()) > 0

async def test_bigsdb_index_references_pubmlst_correctly():
    async with BIGSdbIndex() as bigsdb_index:
        assert (await bigsdb_index.get_bigsdb_api_from_seqdefdb("pubmlst_hinfluenzae_seqdef")) == "https://rest.pubmlst.org"

async def test_bigsdb_index_references_institutpasteur_correctly():
    async with BIGSdbIndex() as bigsdb_index:
        assert (await bigsdb_index.get_bigsdb_api_from_seqdefdb("pubmlst_bordetella_seqdef")) == "https://bigsdb.pasteur.fr/api"


async def test_bigsdb_index_instantiates_correct_profiler():
    sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    async with BIGSdbIndex() as bigsdb_index:
        async with await bigsdb_index.build_profiler_from_seqdefdb("pubmlst_bordetella_seqdef", 3) as profiler:
            profile = await profiler.profile_string(sequence)
            assert profile.clonal_complex == "ST-2 complex"
            assert profile.sequence_type == "1"

async def test_bigsdb_profile_multiple_strings_same_string_twice():
    sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    dummy_sequences = [NamedString("seq1", sequence), NamedString("seq2", sequence)]
    async def generate_async_iterable_sequences():
        for dummy_sequence in dummy_sequences:
            yield [dummy_sequence]
    async with OnlineBIGSdbMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        async for named_profile in dummy_profiler.profile_multiple_strings(generate_async_iterable_sequences()):
            name, profile = named_profile.name, named_profile.mlst_profile
            assert profile is not None
            assert isinstance(profile, MLSTProfile)
            assert profile.clonal_complex == "ST-2 complex"
            assert profile.sequence_type == "1"

async def test_bigsdb_profile_multiple_strings_exactmatch_fail_second_no_stop():
    valid_seq = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    dummy_sequences = [NamedString("seq1", valid_seq), NamedString("should_fail", gene_scrambler(valid_seq, 0.3)), NamedString("seq3", valid_seq)]
    async def generate_async_iterable_sequences():
        for dummy_sequence in dummy_sequences:
            yield [dummy_sequence]
    async with OnlineBIGSdbMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        async for name_profile in dummy_profiler.profile_multiple_strings(generate_async_iterable_sequences(), True):
            name, profile = name_profile.name, name_profile.mlst_profile

            assert profile is not None
            assert isinstance(profile, MLSTProfile)
            if name == "should_fail":
                assert profile.clonal_complex == "unknown"
                assert profile.sequence_type == "unknown"
            else:
                assert profile.clonal_complex == "ST-2 complex"
                assert profile.sequence_type == "1"

async def test_bigsdb_profile_multiple_strings_nonexact_second_no_stop():
    valid_seq = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    dummy_sequences = [NamedString("seq1", valid_seq), NamedString("should_fail", gene_scrambler(valid_seq, 0.3)), NamedString("seq3", valid_seq)]
    async def generate_async_iterable_sequences():
        for dummy_sequence in dummy_sequences:
            yield [dummy_sequence]
    async with OnlineBIGSdbMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        async for named_profile in dummy_profiler.profile_multiple_strings(generate_async_iterable_sequences(), False):
            name, profile = named_profile.name, named_profile.mlst_profile
            if name == "should_fail":
                assert profile is not None
                assert profile.clonal_complex == "unknown"
                assert profile.sequence_type == "unknown"
                assert len(profile.alleles) > 0
            else:
                assert profile is not None
                assert isinstance(profile, MLSTProfile)
                assert profile.clonal_complex == "ST-2 complex"
                assert profile.sequence_type == "1"

async def test_bigsdb_index_get_schemas_for_bordetella():
    async with BIGSdbIndex() as index:
        schemas = await index.get_schemas_for_seqdefdb(seqdef_db_name="pubmlst_bordetella_seqdef")
        assert len(schemas.keys()) > 0
        assert "MLST" in schemas
        assert isinstance(schemas["MLST"], int)

async def test_bigsdb_index_get_databases_has_only_seqdef():
    async with BIGSdbIndex() as index:
        databases = await index.get_known_seqdef_dbs()
        assert len(databases.keys()) > 0
        for database_name in databases.keys():
            assert database_name.endswith("seqdef")
        assert databases["pubmlst_bordetella_seqdef"] == "https://bigsdb.pasteur.fr/api"