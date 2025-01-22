from collections import defaultdict
from contextlib import AbstractAsyncContextManager
from numbers import Number
from typing import Any, AsyncGenerator, AsyncIterable, Collection, Generator, Iterable, Mapping, Sequence, Union

from aiohttp import ClientSession, ClientTimeout

from autobigs.engine.data.structures.genomics import NamedString
from autobigs.engine.data.structures.mlst import Allele, PartialAllelicMatchProfile, MLSTProfile
from autobigs.engine.exceptions.database import NoBIGSdbExactMatchesException, NoBIGSdbMatchesException, NoSuchBIGSdbDatabaseException

class BIGSdbMLSTProfiler(AbstractAsyncContextManager):

    def __init__(self, database_api: str, database_name: str, schema_id: int):
        self._database_name = database_name
        self._schema_id = schema_id
        self._base_url = f"{database_api}/db/{self._database_name}/schemes/{self._schema_id}/"
        self._http_client = ClientSession(self._base_url, timeout=ClientTimeout(10000))

    async def __aenter__(self):
        return self

    async def fetch_mlst_allele_variants(self, sequence_string: str, exact: bool) -> AsyncGenerator[Allele, Any]:
        # See https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes
        uri_path = "sequence"
        response = await self._http_client.post(uri_path, json={
            "sequence": sequence_string,
            "partial_matches": not exact
        })
        sequence_response: dict = await response.json()

        if "exact_matches" in sequence_response:
            # loci -> list of alleles with id and loci
            exact_matches: dict[str, Sequence[dict[str, str]]] = sequence_response["exact_matches"]  
            for allele_loci, alleles in exact_matches.items():
                for allele in alleles:
                    alelle_id = allele["allele_id"]
                    yield Allele(allele_loci=allele_loci, allele_variant=alelle_id, partial_match_profile=None)
        elif "partial_matches" in sequence_response:
            if exact:
                raise NoBIGSdbExactMatchesException(self._database_name, self._schema_id)
            partial_matches: dict[str, dict[str, Union[str, float, int]]] = sequence_response["partial_matches"] 
            for allele_loci, partial_match in partial_matches.items():
                if len(partial_match) <= 0:
                    continue
                partial_match_profile = PartialAllelicMatchProfile(
                    percent_identity=float(partial_match["identity"]),
                    mismatches=int(partial_match["mismatches"]),
                    bitscore=float(partial_match["bitscore"]),
                    gaps=int(partial_match["gaps"])
                )
                yield Allele(
                    allele_loci=allele_loci,
                    allele_variant=str(partial_match["allele"]),
                    partial_match_profile=partial_match_profile
                )
        else:
            raise NoBIGSdbMatchesException(self._database_name, self._schema_id)



    async def fetch_mlst_st(self, alleles: Union[AsyncIterable[Allele], Iterable[Allele]]) -> MLSTProfile:
        uri_path = "designations"
        allele_request_dict: dict[str, list[dict[str, str]]] = defaultdict(list)
        if isinstance(alleles, AsyncIterable):
            async for allele in alleles:
                allele_request_dict[allele.allele_loci].append({"allele": str(allele.allele_variant)})
        else:
            for allele in alleles:
                allele_request_dict[allele.allele_loci].append({"allele": str(allele.allele_variant)})
        request_json = {
            "designations": allele_request_dict
        }
        async with self._http_client.post(uri_path, json=request_json) as response:
            response_json: dict = await response.json()
            allele_map: dict[str, list[Allele]] = defaultdict(list)
            response_json.setdefault("fields", dict())
            schema_fields_returned: dict[str, str] = response_json["fields"]
            schema_fields_returned.setdefault("ST", "unknown")
            schema_fields_returned.setdefault("clonal_complex", "unknown")
            schema_exact_matches: dict = response_json["exact_matches"]
            for exact_match_loci, exact_match_alleles in schema_exact_matches.items():
                for exact_match_allele in exact_match_alleles:
                    allele_map[exact_match_loci].append(Allele(exact_match_loci, exact_match_allele["allele_id"], None))
            if len(allele_map) == 0:
                raise ValueError("Passed in no alleles.")
            return MLSTProfile(dict(allele_map), schema_fields_returned["ST"], schema_fields_returned["clonal_complex"])

    async def profile_string(self, string: str, exact: bool = False) -> MLSTProfile:
        alleles = self.fetch_mlst_allele_variants(string, exact)
        return await self.fetch_mlst_st(alleles)


    async def profile_multiple_strings(self, namedStrings: AsyncIterable[NamedString], exact: bool = False, stop_on_fail: bool = False) -> AsyncGenerator[tuple[str, Union[MLSTProfile, None]], Any]:
        async for named_string in namedStrings:
            try:
                yield (named_string.name, await self.profile_string(named_string.sequence, exact))
            except NoBIGSdbMatchesException as e:
                if stop_on_fail:
                    raise e
                yield (named_string.name, None)

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

    async def build_profiler_from_seqdefdb(self, dbseqdef_name: str, schema_id: int) -> BIGSdbMLSTProfiler:
        return BIGSdbMLSTProfiler(await self.get_bigsdb_api_from_seqdefdb(dbseqdef_name), dbseqdef_name, schema_id)

    async def close(self):
        await self._http_client.close()

    async def __aexit__(self, exc_type, exc_value, traceback):
        await self.close()
    
