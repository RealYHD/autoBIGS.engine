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
    def determine_mlst_allele_variants(self, query_sequence_strings: Union[Iterable[Union[NamedString, str]], Union[NamedString, str]]) -> AsyncGenerator[Union[Allele, tuple[str, Allele]], Any]:
        pass

    @abstractmethod
    async def determine_mlst_st(self, alleles: Union[AsyncIterable[Union[Allele, tuple[str, Allele]]], Iterable[Union[Allele, tuple[str, Allele]]]]) -> Union[MLSTProfile, NamedMLSTProfile]:
        pass

    @abstractmethod
    async def profile_string(self, query_sequence_strings: Iterable[Union[NamedString, str]]) -> Union[NamedMLSTProfile, MLSTProfile]:
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

    async def determine_mlst_allele_variants(self, query_sequence_strings: Union[Iterable[Union[NamedString, str]], Union[NamedString, str]]) -> AsyncGenerator[Union[Allele, tuple[str, Allele]], Any]:
        # See https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes
        uri_path = "sequence"
        if isinstance(query_sequence_strings, str) or isinstance(query_sequence_strings, NamedString):
            query_sequence_strings = [query_sequence_strings]
        for sequence_string in query_sequence_strings:
            async with self._http_client.post(uri_path, json={
                "sequence": sequence_string if isinstance(sequence_string, str) else sequence_string.sequence,
                "partial_matches": True
            }) as response:
                sequence_response: dict = await response.json()

                if "exact_matches" in sequence_response:
                    # loci -> list of alleles with id and loci
                    exact_matches: dict[str, Sequence[dict[str, str]]] = sequence_response["exact_matches"]  
                    for allele_loci, alleles in exact_matches.items():
                        for allele in alleles:
                            alelle_id = allele["allele_id"]
                            result_allele = Allele(allele_locus=allele_loci, allele_variant=alelle_id, partial_match_profile=None)
                            yield result_allele if isinstance(sequence_string, str) else (sequence_string.name, result_allele)
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
                        result_allele = Allele(
                            allele_locus=allele_loci,
                            allele_variant=str(partial_match["allele"]),
                            partial_match_profile=partial_match_profile
                        )
                        yield result_allele if isinstance(sequence_string, str) else (sequence_string.name, result_allele)
                else:
                    raise NoBIGSdbMatchesException(self._database_name, self._schema_id, sequence_string.name if isinstance(sequence_string, NamedString) else None)

    async def determine_mlst_st(self, alleles: Union[AsyncIterable[Union[Allele, tuple[str, Allele]]], Iterable[Union[Allele, tuple[str, Allele]]]]) -> Union[MLSTProfile, NamedMLSTProfile]:
        uri_path = "designations"
        allele_request_dict: dict[str, list[dict[str, str]]] = defaultdict(list)
        names_list = []
        def insert_allele_to_request_dict(allele: Union[Allele, tuple[str, Allele]]):
            if isinstance(allele, Allele):
                allele_val = allele
            else:
                allele_val = allele[1]
                names_list.append(allele[0])
            allele_request_dict[allele_val.allele_locus].append({"allele": str(allele_val.allele_variant)})

        if isinstance(alleles, AsyncIterable):
            async for allele in alleles:
                insert_allele_to_request_dict(allele)
        else:
            for allele in alleles:
                insert_allele_to_request_dict(allele)
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
                allele_set.add(Allele(exact_match_locus, exact_match_alleles[0]["allele_id"], None))
            if len(allele_set) == 0:
                raise ValueError("Passed in no alleles.")
            result_mlst_profile = MLSTProfile(allele_set, schema_fields_returned["ST"], schema_fields_returned["clonal_complex"])
            if len(names_list) > 0:
                result_mlst_profile = NamedMLSTProfile(str(tuple(names_list)), result_mlst_profile)
            return result_mlst_profile

    async def profile_string(self, query_sequence_strings: Iterable[Union[NamedString, str]]) -> Union[NamedMLSTProfile, MLSTProfile]:
        alleles = self.determine_mlst_allele_variants(query_sequence_strings)
        return await self.determine_mlst_st(alleles)

    async def profile_multiple_strings(self, query_named_string_groups: AsyncIterable[Iterable[NamedString]], stop_on_fail: bool = False) -> AsyncGenerator[NamedMLSTProfile, Any]:
        tasks = []
        async for named_strings in query_named_string_groups:
            tasks.append(self.profile_string(named_strings))
            for task in asyncio.as_completed(tasks):
                try:
                    yield await task
                except NoBIGSdbMatchesException as e:
                    if stop_on_fail:
                        raise e
                    causal_name = e.get_causal_query_name()
                    if causal_name is None:
                        raise ValueError("Missing query name despite requiring names.")
                    else:
                        yield NamedMLSTProfile(causal_name, None)

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