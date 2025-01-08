import asyncio
from io import TextIOWrapper
from typing import Any, AsyncGenerator, Generator, Sequence, Union
from Bio import SeqIO

from automlst.engine.data.genomics import NamedString

async def read_fasta(handle: Union[str, TextIOWrapper]) -> AsyncGenerator[NamedString, Any]:
    fasta_sequences = asyncio.to_thread(SeqIO.parse, handle=handle, format="fasta")
    for fasta_sequence in await fasta_sequences:
        yield NamedString(fasta_sequence.id, str(fasta_sequence.seq))