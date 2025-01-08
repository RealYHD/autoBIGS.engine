
import asyncio
import datetime
from typing import Any, AsyncGenerator, AsyncIterable, Iterable, Sequence, Union
from automlst.cli import program
from automlst.engine.data.genomics import NamedString
from automlst.engine.data.mlst import MLSTProfile
from automlst.engine.local.abif import read_abif, reference_consensus_assembly
from automlst.engine.local.csv import write_mlst_profiles_as_csv
from automlst.engine.local.fasta import read_fasta, read_multiple_fastas
from automlst.engine.remote.databases.bigsdb import BIGSdbIndex, BigSDBMLSTProfiler


parser = program.subparsers.add_parser(__name__)

parser.add_argument(
    "--fasta", "-fa", "-fst",
    nargs="+",
    action='extend',
    dest="fastas",
    required=False,
    default=[],
    type=str,
    help="The FASTA files to process. Multiple can be listed."
)

parser.add_argument(
    "seqdefdb",
    help="The BIGSdb seqdef database to use for typing."
)

parser.add_argument(
    "schema",
    type=int,
    help="The BIGSdb seqdef database schema ID (integer) to use for typing."
)

parser.add_argument(
    "out",
    default=f'./{datetime.datetime.now().strftime(r"%Y%m%d%H%M%S")}',
    help="The output CSV name (.csv will be appended)."
)


async def run(args):
    async with BIGSdbIndex() as bigsdb_index:
        gen_strings = read_multiple_fastas(args.fastas)
        async with await bigsdb_index.build_profiler_from_seqdefdb(args.seqdefdb, args.schema) as mlst_profiler:
            mlst_profiles = mlst_profiler.profile_multiple_strings(gen_strings)
            await write_mlst_profiles_as_csv(mlst_profiles, args.out)

def run_asynchronously(args):
    asyncio.run(run(args))

parser.set_defaults(func=run_asynchronously)
