from Bio import SeqIO
from Bio.Align import PairwiseAligner
from pytest import mark
from pytest import fixture
from autobigs.engine.analysis.aligners import AsyncBiopythonPairwiseAlignmentEngine
from autobigs.engine.structures.alignment import PairwiseAlignment

@fixture
def tohamaI_bpertussis_adk():
    return str(SeqIO.read("tests/resources/tohama_I_bpertussis_adk.fasta", format="fasta").seq)

@fixture
def tohamaI_bpertussis_genome():
    return str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", format="fasta").seq)

@fixture
def fdaargos_1560_hinfluenza_adk():
    return str(SeqIO.read("tests/resources/fdaargos_1560_hinfluenza_adk.fasta", format="fasta").seq)

@fixture
def fdaargos_1560_hinfluenza_genome():
    return str(SeqIO.read("tests/resources/fdaargos_1560_hinfluenza.fasta", format="fasta").seq)


@fixture(params=[1, 2])
def dummy_engine(request):
    aligner = PairwiseAligner("blastn")
    aligner.mode = "local"
    with AsyncBiopythonPairwiseAlignmentEngine(aligner, request.param) as engine:
        yield engine

class TestAsyncPairwiseAlignmentEngine:
    async def test_single_alignment_no_errors_single_alignment(self, tohamaI_bpertussis_genome, tohamaI_bpertussis_adk: str, dummy_engine: AsyncBiopythonPairwiseAlignmentEngine):
        dummy_engine.align(tohamaI_bpertussis_genome, tohamaI_bpertussis_adk)
        async for alignment, additional_information in dummy_engine:
            assert isinstance(alignment, PairwiseAlignment)

    async def test_single_alignment_no_errors_multiple(self, tohamaI_bpertussis_genome, tohamaI_bpertussis_adk, fdaargos_1560_hinfluenza_genome, fdaargos_1560_hinfluenza_adk, dummy_engine: AsyncBiopythonPairwiseAlignmentEngine):
        dummy_engine.align(tohamaI_bpertussis_genome, tohamaI_bpertussis_adk)
        dummy_engine.align(fdaargos_1560_hinfluenza_genome, fdaargos_1560_hinfluenza_adk)
        async for alignment, additional_information in dummy_engine:
            assert isinstance(alignment, PairwiseAlignment)