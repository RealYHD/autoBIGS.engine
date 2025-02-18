from abc import abstractmethod
import asyncio
from collections import defaultdict
from contextlib import AbstractAsyncContextManager
import csv
from os import path
import os
import shutil
import tempfile
from typing import Any, AsyncGenerator, AsyncIterable, Iterable, Mapping, Sequence, Set, Union

from aiohttp import ClientSession, ClientTimeout

from autobigs.engine.reading import read_fasta
from autobigs.engine.structures.alignment import PairwiseAlignment
from autobigs.engine.structures.genomics import NamedString
from autobigs.engine.structures.mlst import Allele, NamedMLSTProfile, AlignmentStats, MLSTProfile
from autobigs.engine.exceptions.database import NoBIGSdbExactMatchesException, NoBIGSdbMatchesException, NoSuchBIGSdbDatabaseException

from Bio.Align import PairwiseAligner

class BIGSdbMLSTProfiler(AbstractAsyncContextManager):

    @abstractmethod
    def determine_mlst_allele_variants(self, query_sequence_strings: Iterable[str]) -> AsyncGenerator[Allele, Any]:
        pass

    @abstractmethod
    async def determine_mlst_st(self, alleles: Union[AsyncIterable[Allele], Iterable[Allele]]) -> MLSTProfile:
        pass

    @abstractmethod
    async def profile_string(self, query_sequence_strings: Iterable[str]) -> MLSTProfile:
        pass

    @abstractmethod
    def profile_multiple_strings(self, query_named_string_groups: AsyncIterable[Iterable[NamedString]], stop_on_fail: bool = False) -> AsyncGenerator[NamedMLSTProfile, Any]:
        pass

    @abstractmethod
    async def close(self):
        pass

class RemoteBIGSdbMLSTProfiler(BIGSdbMLSTProfiler):

    def __init__(self, database_api: str, database_name: str, schema_id: int):
        self._database_name = database_name
        self._schema_id = schema_id
        self._base_url = f"{database_api}/db/{self._database_name}/schemes/{self._schema_id}/"
        self._http_client = ClientSession(self._base_url, timeout=ClientTimeout(60))

    async def __aenter__(self):
        return self

    async def determine_mlst_allele_variants(self, query_sequence_strings: Union[Iterable[str], str]) -> AsyncGenerator[Allele, Any]:
        # See https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes
        uri_path = "sequence"
        if isinstance(query_sequence_strings, str):
            query_sequence_strings = [query_sequence_strings]
        for sequence_string in query_sequence_strings:
            async with self._http_client.post(uri_path, json={
                "sequence": sequence_string,
                "partial_matches": True
            }) as response:
                sequence_response: dict = await response.json()

                if "exact_matches" in sequence_response:
                    # loci -> list of alleles with id and loci
                    exact_matches: dict[str, Sequence[dict[str, str]]] = sequence_response["exact_matches"]  
                    for allele_loci, alleles in exact_matches.items():
                        for allele in alleles:
                            alelle_id = allele["allele_id"]
                            yield Allele(allele_locus=allele_loci, allele_variant=alelle_id, partial_match_profile=None)
                elif "partial_matches" in sequence_response:
                    partial_matches: dict[str, dict[str, Union[str, float, int]]] = sequence_response["partial_matches"] 
                    for allele_loci, partial_match in partial_matches.items():
                        if len(partial_match) <= 0:
                            continue
                        partial_match_profile = AlignmentStats(
                            percent_identity=float(partial_match["identity"]),
                            mismatches=int(partial_match["mismatches"]),
                            gaps=int(partial_match["gaps"]),
                            match_metric=int(partial_match["bitscore"])
                        )
                        yield Allele(
                            allele_locus=allele_loci,
                            allele_variant=str(partial_match["allele"]),
                            partial_match_profile=partial_match_profile
                        )
                else:
                    raise NoBIGSdbMatchesException(self._database_name, self._schema_id)

    async def determine_mlst_st(self, alleles: Union[AsyncIterable[Allele], Iterable[Allele]]) -> MLSTProfile:
        uri_path = "designations"
        allele_request_dict: dict[str, list[dict[str, str]]] = defaultdict(list)
        if isinstance(alleles, AsyncIterable):
            async for allele in alleles:
                allele_request_dict[allele.allele_locus].append({"allele": str(allele.allele_variant)})
        else:
            for allele in alleles:
                allele_request_dict[allele.allele_locus].append({"allele": str(allele.allele_variant)})
        request_json = {
            "designations": allele_request_dict
        }
        async with self._http_client.post(uri_path, json=request_json) as response:
            response_json: dict = await response.json()
            allele_set: Set[Allele] = set()
            response_json.setdefault("fields", dict())
            schema_fields_returned: dict[str, str] = response_json["fields"]
            schema_fields_returned.setdefault("ST", "unknown")
            schema_fields_returned.setdefault("clonal_complex", "unknown")
            schema_exact_matches: dict = response_json["exact_matches"]
            for exact_match_locus, exact_match_alleles in schema_exact_matches.items():
                if len(exact_match_alleles) > 1:
                    raise ValueError(f"Unexpected number of alleles returned for exact match (Expected 1, retrieved {len(exact_match_alleles)})")
                allele_set.add(Allele(exact_match_locus, exact_match_alleles[0]["allele_id"], None))
            if len(allele_set) == 0:
                raise ValueError("Passed in no alleles.")
            return MLSTProfile(allele_set, schema_fields_returned["ST"], schema_fields_returned["clonal_complex"])

    async def profile_string(self, query_sequence_strings: Iterable[str]) -> MLSTProfile:
        alleles = self.determine_mlst_allele_variants(query_sequence_strings)
        return await self.determine_mlst_st(alleles)

    async def profile_multiple_strings(self, query_named_string_groups: AsyncIterable[Iterable[NamedString]], stop_on_fail: bool = False) -> AsyncGenerator[NamedMLSTProfile, Any]:
        async for named_strings in query_named_string_groups:
            names: list[str] = list()
            sequences: list[str] = list()
            for named_string in named_strings:
                names.append(named_string.name)
                sequences.append(named_string.sequence)
            try:
                yield NamedMLSTProfile("-".join(names), (await self.profile_string(sequences)))
            except NoBIGSdbMatchesException as e:
                if stop_on_fail:
                    raise e
                yield NamedMLSTProfile("-".join(names), None)

    async def close(self):
        await self._http_client.close()

    async def __aexit__(self, exc_type, exc_value, traceback):
        await self.close()

class BIGSdbIndex(AbstractAsyncContextManager):
    KNOWN_BIGSDB_APIS = {
        "https://bigsdb.pasteur.fr/api",
        "https://rest.pubmlst.org"
    }

    def __init__(self):
        self._http_client = ClientSession()
        self._known_seqdef_dbs_origin: Union[Mapping[str, str], None] = None
        self._seqdefdb_schemas: dict[str, Union[Mapping[str, int], None]] = dict()
        super().__init__()

    async def __aenter__(self):
        return self
    
    async def get_known_seqdef_dbs(self, force: bool = False) -> Mapping[str, str]:
        if self._known_seqdef_dbs_origin is not None and not force:
            return self._known_seqdef_dbs_origin
        known_seqdef_dbs = dict()
        for known_bigsdb in BIGSdbIndex.KNOWN_BIGSDB_APIS:
            async with self._http_client.get(f"{known_bigsdb}/db") as response:
                response_json_databases = await response.json()
                for database_group in response_json_databases:
                    for database_info in database_group["databases"]:
                        if str(database_info["name"]).endswith("seqdef"):
                            known_seqdef_dbs[database_info["name"]] = known_bigsdb
        self._known_seqdef_dbs_origin = dict(known_seqdef_dbs)
        return self._known_seqdef_dbs_origin

    async def get_bigsdb_api_from_seqdefdb(self, seqdef_db_name: str) -> str:
        known_databases = await self.get_known_seqdef_dbs()
        if seqdef_db_name not in known_databases:
            raise NoSuchBIGSdbDatabaseException(seqdef_db_name)
        return known_databases[seqdef_db_name]     

    async def get_schemas_for_seqdefdb(self, seqdef_db_name: str, force: bool = False) -> Mapping[str, int]:
        if seqdef_db_name in self._seqdefdb_schemas and not force:
            return self._seqdefdb_schemas[seqdef_db_name] # type: ignore since it's guaranteed to not be none by conditional
        uri_path = f"{await self.get_bigsdb_api_from_seqdefdb(seqdef_db_name)}/db/{seqdef_db_name}/schemes"
        async with self._http_client.get(uri_path) as response: 
            response_json = await response.json()
            schema_descriptions: Mapping[str, int] = dict()
            for scheme_definition in response_json["schemes"]:
                scheme_id: int = int(str(scheme_definition["scheme"]).split("/")[-1])
                scheme_desc: str = scheme_definition["description"]
                schema_descriptions[scheme_desc] = scheme_id
            self._seqdefdb_schemas[seqdef_db_name] = schema_descriptions
            return self._seqdefdb_schemas[seqdef_db_name] # type: ignore

    async def build_profiler_from_seqdefdb(self, local: bool, dbseqdef_name: str, schema_id: int) -> BIGSdbMLSTProfiler:
        return get_BIGSdb_MLST_profiler(local, await self.get_bigsdb_api_from_seqdefdb(dbseqdef_name), dbseqdef_name, schema_id)

    async def close(self):
        await self._http_client.close()

    async def __aexit__(self, exc_type, exc_value, traceback):
        await self.close()

def get_BIGSdb_MLST_profiler(local: bool, database_api: str, database_name: str, schema_id: int):
    if local:
        raise NotImplementedError()
    return RemoteBIGSdbMLSTProfiler(database_api=database_api, database_name=database_name, schema_id=schema_id)