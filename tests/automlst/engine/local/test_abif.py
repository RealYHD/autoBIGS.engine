import os

from automlst.engine.data.local.abif import read_abif

async def test_load_sanger_sequence_has_data():
    assert os.path.exists("tests/resources/1I1_F_P1815443_047.ab1")
    result_data = await read_abif("tests/resources/1I1_F_P1815443_047.ab1")
    assert result_data is not None
