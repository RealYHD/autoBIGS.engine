import asyncio
from io import TextIOWrapper
from typing import Any, AsyncGenerator, Iterable, Union
from Bio import SeqIO

from autobigs.engine.structures.genomics import NamedString

async def read_fasta(handle: Union[str, TextIOWrapper]) -> AsyncGenerator[NamedString, Any]:
    fasta_sequences = asyncio.to_thread(SeqIO.parse, handle=handle, format="fasta")
    for fasta_sequence in await fasta_sequences:
        yield NamedString(fasta_sequence.id, str(fasta_sequence.seq))

async def read_multiple_fastas(handles: Iterable[Union[str, TextIOWrapper]]) -> AsyncGenerator[NamedString, Any]:
    for handle in handles:
        async for named_seq in read_fasta(handle):
            yield named_seq