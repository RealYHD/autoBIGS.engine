from collections import defaultdict
from contextlib import AbstractAsyncContextManager
import re
from typing import Any, AsyncGenerator, AsyncIterable, Generator, Iterable, Sequence, Union
from aiohttp import ClientSession, ClientTimeout
from automlst.engine.data.MLST import Allele, MLSTProfile
from automlst.engine.data.genomics import NamedString

class InstitutPasteurProfiler(AbstractAsyncContextManager):

    async def __aenter__(self):
        return self


    def __init__(self, database_name: str):
        self._base_url = f"https://bigsdb.pasteur.fr/api/db/{database_name}/"
        self._http_client = ClientSession(self._base_url, timeout=ClientTimeout(10000))

    async def fetch_mlst_allele_variants(self, sequence_string: str) -> AsyncGenerator[Allele, Any]:
        # See https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes
        uri_path = "schemes/3/sequence"
        response = await self._http_client.post(uri_path, json={
            "sequence": sequence_string
        })
        sequence_response: dict = await response.json()
        exact_matches: dict[str, Sequence[dict[str, str]]] = sequence_response["exact_matches"]  
        for allele_loci, alleles in exact_matches.items():
            for allele in alleles:
                alelle_id = allele["allele_id"]
                yield Allele(allele_loci=allele_loci, allele_variant=alelle_id)

    async def fetch_mlst_st(self, alleles: Union[AsyncIterable[Allele], Iterable[Allele]]) -> MLSTProfile:
        uri_path = "schemes/3/designations"
        allele_request_dict: dict[str, list[dict[str, str]]] = defaultdict(list)
        if isinstance(alleles, AsyncIterable):
            async for allele in alleles:
                allele_request_dict[allele.allele_loci].append({"allele": str(allele.allele_variant)})
        else:
            for allele in alleles:
                allele_request_dict[allele.allele_loci].append({"allele": str(allele.allele_variant)})
        response = await self._http_client.post(uri_path, json={
            "designations": allele_request_dict
        })
        response_json = await response.json()
        schema_fields_returned = response_json["fields"]
        schema_exact_matches = response_json["exact_matches"]
        allele_map: dict[str, list[Allele]] = defaultdict(list)
        for exact_match_loci, exact_match_alleles in schema_exact_matches.items():
            for exact_match_allele in exact_match_alleles:
                allele_map[exact_match_loci].append(Allele(exact_match_loci, exact_match_allele["allele_id"]))
        return MLSTProfile(allele_map, schema_fields_returned["ST"], schema_fields_returned["clonal_complex"])

    async def profile_string(self, string: str) -> MLSTProfile:
        alleles = self.fetch_mlst_allele_variants(string)
        return await self.fetch_mlst_st(alleles)


    async def close(self):
        await self._http_client.close()

    async def __aexit__(self, exc_type, exc_value, traceback):
        await self.close()