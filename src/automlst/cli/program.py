import argparse
import asyncio
import datetime
from os import path
import os

from automlst.engine.data.genomics import NamedString
from automlst.engine.local.abif import read_abif
from automlst.engine.local.csv import write_mlst_profiles_as_csv
from automlst.engine.local.fasta import read_fasta
from automlst.engine.remote.databases.bigsdb import BIGSdbIndex

root_parser = argparse.ArgumentParser()
subparsers = root_parser.add_subparsers(required=True)

def run():
    args = root_parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    run()