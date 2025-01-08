from Bio import SeqIO
from automlst.engine.data.mlst import Allele, MLSTProfile
from automlst.engine.remote.databases.institutpasteur.mlst import InstitutPasteurProfiler


async def test_profiling_results_in_exact_matches_when_exact():
    sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    async with InstitutPasteurProfiler(database_name="pubmlst_bordetella_seqdef") as dummy_profiler:
        exact_matches = dummy_profiler.fetch_mlst_allele_variants(schema_id=3, sequence_string=sequence)
        targets_left = {"adk", "fumC", "glyA", "tyrB", "icd", "pepA", "pgm"}
        async for exact_match in exact_matches:
            assert isinstance(exact_match, Allele)
            assert exact_match.allele_variant == '1' # All of Tohama I has allele id I
            targets_left.remove(exact_match.allele_loci)

        assert len(targets_left) == 0

async def test_profiling_results_in_correct_st():
    dummy_alleles = [
        Allele("adk", "1"),
        Allele("fumC", "1"),
        Allele("glyA", "1"),
        Allele("tyrB", "1"),
        Allele("icd", "1"),
        Allele("pepA", "1"),
        Allele("pgm", "1"),
    ]
    async with InstitutPasteurProfiler(database_name="pubmlst_bordetella_seqdef") as dummy_profiler:
        mlst_st_data = await dummy_profiler.fetch_mlst_st(3, dummy_alleles)
        assert mlst_st_data is not None
        assert isinstance(mlst_st_data, MLSTProfile)
        assert mlst_st_data.clonal_complex == "ST-2 complex"
        assert mlst_st_data.sequence_type == "1"

async def test_sequence_profiling_is_correct():
    sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    dummy_alleles = [
        Allele("adk", "1"),
        Allele("fumC", "1"),
        Allele("glyA", "1"),
        Allele("tyrB", "1"),
        Allele("icd", "1"),
        Allele("pepA", "1"),
        Allele("pgm", "1"),
    ]
    async with InstitutPasteurProfiler(database_name="pubmlst_bordetella_seqdef") as dummy_profiler:
        profile = await dummy_profiler.profile_string(3, sequence)
        assert profile is not None
        assert isinstance(profile, MLSTProfile)
        assert profile.clonal_complex == "ST-2 complex"
        assert profile.sequence_type == "1"