import argparse
import asyncio
import datetime
from os import path
import os

from automlst.cli import aggregated
from automlst.engine.data.genomics import NamedString
from automlst.engine.local.abif import read_abif
from automlst.engine.local.csv import write_mlst_profiles_as_csv
from automlst.engine.local.fasta import read_fasta


parser = argparse.ArgumentParser()
parser.add_argument(
    "--run-name", "-name",
    dest="run_name",
    required=False,
    default=datetime.datetime.now().strftime(r"%Y%m%d%H%M%S"),
    type=str,
    help="The name of the run. Will use a date and time string if not provided."
)
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
    "--abif", "-abi", "-ab1",
    action='extend',
    dest="abifs",
    required=False,
    default=[],
    type=str,
    help="The ABIF files to process. Multiple can be listed."
)
parser.add_argument(
    "--ncbi-assembly-reference", "-refncbi",
    dest="ncbi_assembly_reference",
    required=False,
    default=None,
    type=str,
    help="The NCBI GenBank accession ID for the consensus assembly. Either this argument, or the path equivalent must be given if ABIF files are used."
)
parser.add_argument(
    "--assembly-reference", "-ref",
    dest="assembly_reference",
    required=False,
    default=None,
    type=str,
    help="The path to the FASTA sequence to be used as a reference for consensus building."
)
parser.add_argument(
    "--institut-pasteur-mlst",
    "-ipdbmlst",
    dest="institut_pasteur_db",
    required=False,
    default=None,
    type=str,
    help="The Institut Pasteur MLST database to use."
)
parser.add_argument(
    "out",
    default="./.",
    help="The output folder. Files will be named by the provided (or default) run name. Outputs will be automatically generated depending on which arguments are used."
)


def run():
    args = parser.parse_args()
    gen_strings = aggregated.aggregate_sequences(args.fastas, args.abifs)
    os.makedirs(args.out, exist_ok=True)
    if args.institut_pasteur_db is not None:
        mlst_profiles = aggregated.profile_all_genetic_strings(
            gen_strings, args.institut_pasteur_db)
        asyncio.run(write_mlst_profiles_as_csv(
            asyncio.run(mlst_profiles), str(path.join(args.out, "MLST_" + args.run_name + ".csv"))))


if __name__ == "__main__":
    run()