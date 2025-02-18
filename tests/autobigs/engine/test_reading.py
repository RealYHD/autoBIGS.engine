from autobigs.engine.reading import read_fasta


async def test_fasta_reader_not_none():
    named_strings = await read_fasta("tests/resources/tohama_I_bpertussis.fasta")
    for named_string in named_strings:
        assert named_string.name == "BX470248.1"
