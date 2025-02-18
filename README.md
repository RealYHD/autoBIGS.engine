# autoBIGS.engine

A python library implementing common BIGSdb MLST schemes and databases accesses for the purpose of typing sequences automatically. Implementation follows the RESTful API outlined by the official [BIGSdb documentation](https://bigsdb.readthedocs.io/en/latest/rest.html) up to `V1.50.0`.


## Features

Briefly, this library can:
- Import multiple `FASTA` files
- Fetch the available BIGSdb databases that is currently live and available
- Fetch the available BIGSdb database schemas for a given MLST database
- Retrieve exact/non-exact MLST allele variant IDs based off a sequence
- Retrieve MLST sequence type IDs based off a sequence
- Output all results to a single CSV

Furthermore, this library is highly asynchronous where any potentially blocking operation, ranging from parsing FASTAs to performing HTTP requests are at least asynchronous, if not fully multithreaded.

## Usage

This library can be installed through pip. Learn how to [setup and install pip first](https://pip.pypa.io/en/stable/installation/).

Then, it's as easy as running `pip install autobigs-engine` in any terminal that has pip in it's path (any terminal where `pip --version` returns a valid version and install path).

### CLI usage

This is a independent python library and thus does not have any form of direct user interface. One way of using it could be to create your own Python script that makes calls to this libraries functions. Alternatively, you may use `autobigs-cli`, a `Python` package that implements a CLI for calling this library.

## Versioning

the autoBIGS project follows [semantic versioning](https://semver.org/) where the three numbers may be interpreted as MAJOR.MINOR.PATCH.

Note regarding major version 0 ([spec item 4](https://semver.org/#spec-item-4)), the following adaptation of semantic versioning definition is as follows:

1. Given x.Y.z, Y is only incremented when a backwards incompatible change is made.

2. Given x.y.Z, Z is only incremented when a backwards compatible change is made.

Versions of autoBIGS items with a major version number of 0 will introduce numerous changes and patches. As such, changes between such versions should be considered highly variable.