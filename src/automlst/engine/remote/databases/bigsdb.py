from collections import defaultdict
from contextlib import AbstractAsyncContextManager
from typing import Any, AsyncGenerator, AsyncIterable, Collection, Generator, Iterable, Mapping, Sequence, Union

from aiohttp import ClientSession, ClientTimeout

from automlst.engine.data.genomics import NamedString
from automlst.engine.data.mlst import Allele, MLSTProfile

class BIGSdbMLSTProfiler(AbstractAsyncContextManager):

    def __init__(self, database_api: str, database_name: str, schema_id: int):
        self._database_name = database_name
        self._schema_id = schema_id
        self._base_url = f"{database_api}/db/{self._database_name}/schemes/{self._schema_id}/"
        self._http_client = ClientSession(self._base_url, timeout=ClientTimeout(10000))

    async def __aenter__(self):
        return self

    async def fetch_mlst_allele_variants(self, sequence_string: str) -> AsyncGenerator[Allele, Any]:
        # See https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes
        uri_path = "sequence"
        response = await self._http_client.post(uri_path, json={
            "sequence": sequence_string
        })
        sequence_response: dict = await response.json()
        if "exact_matches" not in sequence_response:
            # TODO throw exception for not finding matches.
            pass

        if "exact_matches" not in sequence_response:
            raise ValueError(f"Unable to find exact matches in \"{self._database_name}\" under schema ID \"{self._schema_id}\".")
        exact_matches: dict[str, Sequence[dict[str, str]]] = sequence_response["exact_matches"]  
        for allele_loci, alleles in exact_matches.items():
            for allele in alleles:
                alelle_id = allele["allele_id"]
                yield Allele(allele_loci=allele_loci, allele_variant=alelle_id)

    async def fetch_mlst_st(self, alleles: AsyncIterable[Allele]) -> MLSTProfile:
        uri_path = "designations"
        allele_request_dict: dict[str, list[dict[str, str]]] = defaultdict(list)
        async for allele in alleles:
            allele_request_dict[allele.allele_loci].append({"allele": str(allele.allele_variant)})

        request_json = {
            "designations": allele_request_dict
        }
        async with self._http_client.post(uri_path, json=request_json) as response:
            response_json = await response.json()
            if "fields" not in response_json:
                # TODO raise exception about invalid parameters or no exact parameterization found
                pass
            schema_fields_returned = response_json["fields"]
            schema_exact_matches: dict = response_json["exact_matches"]
            allele_map: dict[str, list[Allele]] = defaultdict(list)
            for exact_match_loci, exact_match_alleles in schema_exact_matches.items():
                for exact_match_allele in exact_match_alleles:
                    allele_map[exact_match_loci].append(Allele(exact_match_loci, exact_match_allele["allele_id"]))
            return MLSTProfile(allele_map, schema_fields_returned["ST"], schema_fields_returned["clonal_complex"])

    async def profile_string(self, string: str) -> MLSTProfile:
        alleles = self.fetch_mlst_allele_variants(string)
        return await self.fetch_mlst_st(alleles)


    async def profile_multiple_strings(self, namedStrings: AsyncIterable[NamedString]) -> AsyncGenerator[tuple[str, MLSTProfile], Any]:
        async for named_string in namedStrings:
            yield (named_string.name, await self.profile_string(named_string.sequence))


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
            raise ValueError(f"The database \"{seqdef_db_name}\" could not be found.")
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
    
