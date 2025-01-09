from argparse import ArgumentParser
import asyncio
from automlst.engine.remote.databases.bigsdb import BIGSdbIndex


def setup_parser(parser: ArgumentParser):
    parser.description = "Fetches the latest BIGSdb MLST database definitions."
    parser.add_argument(
        "--retrieve-bigsdbs", "-l",
        action="store_true",
        dest="list_dbs",
        required=False,
        default=False,
        help="Lists all known BIGSdb MLST databases (fetched from known APIs and cached)."
    )

    parser.add_argument(
        "--retrieve-bigsdb-schemas", "-lschemas",
        nargs="+",
        action="extend",
        dest="list_bigsdb_schemas",
        required=False,
        default=[],
        type=str,
        help="Lists the known schema IDs for a given BIGSdb sequence definition database name. The name, and then the ID of the schema is given."
    )

    parser.set_defaults(func=run_asynchronously)

async def run(args):
    async with BIGSdbIndex() as bigsdb_index:
        if args.list_dbs:
            known_seqdef_dbs = await bigsdb_index.get_known_seqdef_dbs(force=False)
            print("\n".join(known_seqdef_dbs.keys()))

        for bigsdb_schema_name in args.list_bigsdb_schemas:
            schemas = await bigsdb_index.get_schemas_for_seqdefdb(bigsdb_schema_name)
            for schema_desc, schema_id in schemas.items():
                print(f"{schema_desc}: {schema_id}")

def run_asynchronously(args):
    asyncio.run(run(args))

