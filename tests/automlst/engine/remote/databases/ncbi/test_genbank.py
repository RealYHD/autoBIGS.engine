from automlst.engine.remote.databases.ncbi.genbank import fetch_ncbi_genbank


async def test_fetch_ncbi_genbank_with_id_works():
    assert len((await fetch_ncbi_genbank("CP011448.1")).sequence) > 0 