# autoMLST

A CLI/library for rapidly performing MLST typing via accessing pubMLST and InstitutPasteur MSLT databases.

# Components

## automlst.cli

The command line interface, sets up very minimal and mostly makes calls to the library. Uses argparse and is split into two parts:

- `automlst info`: Provides user information on available databases to pull from, and the schemas available.
- `automlst exactmatch`: Provides users the ability to request exact match results from a given database and schema
