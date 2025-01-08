from os import path
from typing import Any, AsyncGenerator, AsyncIterable, Iterable, Sequence
from automlst.engine.data.mlst import MLSTProfile
from automlst.engine.data.genomics import NamedString
from automlst.engine.local.abif import read_abif
from automlst.engine.local.fasta import read_fasta
from automlst.engine.remote.databases.institutpasteur.mlst import InstitutPasteurProfiler


async def aggregate_sequences(fastas: Iterable[str], abifs: Iterable[str]) -> AsyncGenerator[str, Any]:
    for fasta_path in fastas:
        async for fasta in read_fasta(fasta_path):
            yield fasta.sequence
    for abif_path in abifs:
        abif_data = await read_abif(abif_path)
        yield "".join(abif_data.sequence)

async def profile_all_genetic_strings(strings: AsyncIterable[str], database_name: str) -> Sequence[MLSTProfile]:
    profiles = list()
    async with InstitutPasteurProfiler(database_name=database_name) as profiler:
        async for string in strings:
            profiles.append(await profiler.profile_string(string))
    return profiles