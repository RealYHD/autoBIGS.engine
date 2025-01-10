import argparse
import asyncio
import datetime
from os import path
import os

from automlst.cli import info, st
from automlst.cli.meta import get_module_base_name
from automlst.engine.data.genomics import NamedString
from automlst.engine.local.abif import read_abif
from automlst.engine.local.csv import write_mlst_profiles_as_csv
from automlst.engine.local.fasta import read_fasta
from automlst.engine.remote.databases.bigsdb import BIGSdbIndex

root_parser = argparse.ArgumentParser()
subparsers = root_parser.add_subparsers(required=True)

info.setup_parser(subparsers.add_parser(get_module_base_name(info.__name__)))
st.setup_parser(subparsers.add_parser(get_module_base_name(st.__name__)))


def run():
    args = root_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    run()