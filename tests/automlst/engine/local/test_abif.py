import os

from automlst.engine.local.abif import read_abif, reference_consensus_assembly

async def test_load_sanger_sequence_has_data():
    assert os.path.exists("tests/resources/1I1_F_P1815443_047.ab1")
    result_data = await read_abif("tests/resources/1I1_F_P1815443_047.ab1")
    assert result_data is not None

async def test_consensus_assembly_with_ncbi():
    consensus = reference_consensus_assembly("ON685494.1", [await read_abif("tests/resources/1I1_F_P1815443_047.ab1"), await read_abif("tests/resources/1I1_R_P1815443_094.ab1")])
    # TODO complete implementing this