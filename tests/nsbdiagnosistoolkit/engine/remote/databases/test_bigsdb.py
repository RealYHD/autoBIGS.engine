from Bio import SeqIO
from automlst.engine.data.mlst import Allele, MLSTProfile
from automlst.engine.remote.databases.bigsdb import BIGSdbIndex, BigSDBMLSTProfiler


async def test_institutpasteur_profiling_results_in_exact_matches_when_exact():
    sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    async with BigSDBMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        exact_matches = dummy_profiler.fetch_mlst_allele_variants(sequence_string=sequence)
        targets_left = {"adk", "fumC", "glyA", "tyrB", "icd", "pepA", "pgm"}
        async for exact_match in exact_matches:
            assert isinstance(exact_match, Allele)
            assert exact_match.allele_variant == '1' # All of Tohama I has allele id I
            targets_left.remove(exact_match.allele_loci)

        assert len(targets_left) == 0

async def test_institutpasteur_profiling_results_in_correct_mlst_st():
    async def dummy_allele_generator():
        dummy_alleles = [
        Allele("adk", "1"),
        Allele("fumC", "1"),
        Allele("glyA", "1"),
        Allele("tyrB", "1"),
        Allele("icd", "1"),
        Allele("pepA", "1"),
        Allele("pgm", "1"),
        ]
        for dummy_allele in dummy_alleles:
            yield dummy_allele
    async with BigSDBMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        mlst_st_data = await dummy_profiler.fetch_mlst_st(dummy_allele_generator())
        assert mlst_st_data is not None
        assert isinstance(mlst_st_data, MLSTProfile)
        assert mlst_st_data.clonal_complex == "ST-2 complex"
        assert mlst_st_data.sequence_type == "1"

async def test_institutpasteur_sequence_profiling_is_correct():
    sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    async with BigSDBMLSTProfiler(database_api="https://bigsdb.pasteur.fr/api", database_name="pubmlst_bordetella_seqdef", schema_id=3) as dummy_profiler:
        profile = await dummy_profiler.profile_string(sequence)
        assert profile is not None
        assert isinstance(profile, MLSTProfile)
        assert profile.clonal_complex == "ST-2 complex"
        assert profile.sequence_type == "1"


async def test_pubmlst_profiling_results_in_exact_matches_when_exact():
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
    async with BigSDBMLSTProfiler(database_api="https://rest.pubmlst.org/", database_name="pubmlst_hinfluenzae_seqdef", schema_id=1) as dummy_profiler:
        exact_matches = dummy_profiler.fetch_mlst_allele_variants(sequence_string=sequence)
        async for exact_match in exact_matches:
            assert isinstance(exact_match, Allele)
            dummy_alleles.remove(exact_match)

        assert len(dummy_alleles) == 0

async def test_pubmlst_profiling_results_in_correct_st():
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
    async with BigSDBMLSTProfiler(database_api="https://rest.pubmlst.org/", database_name="pubmlst_hinfluenzae_seqdef", schema_id=1) as dummy_profiler:
        mlst_st_data = await dummy_profiler.fetch_mlst_st(generate_dummy_targets())
        assert mlst_st_data is not None
        assert isinstance(mlst_st_data, MLSTProfile)
        assert mlst_st_data.clonal_complex == "ST-3 complex"
        assert mlst_st_data.sequence_type == "3"

async def test_pubmlst_sequence_profiling_is_correct():
    sequence = str(SeqIO.read("tests/resources/FDAARGOS_1560.fasta", "fasta").seq)
    async with BigSDBMLSTProfiler(database_api="https://rest.pubmlst.org/", database_name="pubmlst_hinfluenzae_seqdef", schema_id=1) as dummy_profiler:
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