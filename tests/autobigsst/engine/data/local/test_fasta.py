from autobigsst.engine.data.local.fasta import read_fasta


async def test_fasta_reader_not_none():
    named_strings = read_fasta("tests/resources/tohama_I_bpertussis.fasta")
    async for named_string in named_strings:
        assert named_string.name == "BX470248.1"
