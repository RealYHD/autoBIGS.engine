import asyncio
from io import TextIOWrapper
from typing import Any, AsyncGenerator, Iterable, Union
from Bio import SeqIO

from autobigs.engine.structures.genomics import NamedString

async def read_fasta(handle: Union[str, TextIOWrapper]) -> Iterable[NamedString]:
    fasta_sequences = asyncio.to_thread(SeqIO.parse, handle=handle, format="fasta")
    results = []
    for fasta_sequence in await fasta_sequences:
        results.append(NamedString(fasta_sequence.id, str(fasta_sequence.seq)))
    return results

async def read_multiple_fastas(handles: Iterable[Union[str, TextIOWrapper]]) -> AsyncGenerator[Iterable[NamedString], Any]:
    for handle in handles:
        yield await read_fasta(handle)