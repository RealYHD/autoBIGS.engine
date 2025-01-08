import asyncio
from Bio import SeqIO
from automlst.engine.data.mlst import Allele, MLSTProfile
from automlst.engine.remote.databases.pubmlst.mlst import PubMLSTProfiler


async def test_profiling_results_in_exact_matches_when_exact():
    dummy_alleles = {
        Allele("adk", "1"),
        Allele("atpG", "1"),
        Allele("frdB", "1"),
        Allele("fucK", "1"),
        Allele("mdh", "1"),
        Allele("pgi", "1"),
        Allele("recA", "5"),
    }
    sequence = str(SeqIO.read("tests/resources/FDAARGOS_1560.fasta", "fasta").seq)
    async with PubMLSTProfiler(database_name="pubmlst_hinfluenzae_seqdef") as dummy_profiler:
        exact_matches = dummy_profiler.fetch_mlst_allele_variants(schema_id=1, sequence_string=sequence)
        async for exact_match in exact_matches:
            assert isinstance(exact_match, Allele)
            dummy_alleles.remove(exact_match)

        assert len(dummy_alleles) == 0

async def test_profiling_results_in_correct_st():
    async def generate_dummy_targets():
        dummy_alleles = [
                Allele("adk", "1"),
                Allele("atpG", "1"),
                Allele("frdB", "1"),
                Allele("fucK", "1"),
                Allele("mdh", "1"),
                Allele("pgi", "1"),
                Allele("recA", "5"),
            ]
        for dummy_allele in dummy_alleles:
            yield dummy_allele
    async with PubMLSTProfiler(database_name="pubmlst_hinfluenzae_seqdef") as dummy_profiler:
        mlst_st_data = await dummy_profiler.fetch_mlst_st(1, generate_dummy_targets())
        assert mlst_st_data is not None
        assert isinstance(mlst_st_data, MLSTProfile)
        assert mlst_st_data.clonal_complex == "ST-3 complex"
        assert mlst_st_data.sequence_type == "3"

async def test_sequence_profiling_is_correct():
    sequence = str(SeqIO.read("tests/resources/FDAARGOS_1560.fasta", "fasta").seq)
    async with PubMLSTProfiler(database_name="pubmlst_hinfluenzae_seqdef") as dummy_profiler:
        profile = await dummy_profiler.profile_string(1, sequence)
        assert profile is not None
        assert isinstance(profile, MLSTProfile)
        assert profile.clonal_complex == "ST-3 complex"
        assert profile.sequence_type == "3"