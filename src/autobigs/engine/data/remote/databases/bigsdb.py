from abc import abstractmethod
from collections import defaultdict
from contextlib import AbstractAsyncContextManager
import csv
from os import path
from typing import Any, AsyncGenerator, AsyncIterable, Iterable, Mapping, Sequence, Union

from aiohttp import ClientSession, ClientTimeout

from autobigs.engine.data.local.fasta import read_fasta
from autobigs.engine.data.structures.genomics import NamedString
from autobigs.engine.data.structures.mlst import Allele, NamedMLSTProfile, PartialAllelicMatchProfile, MLSTProfile
from autobigs.engine.exceptions.database import NoBIGSdbExactMatchesException, NoBIGSdbMatchesException, NoSuchBIGSdbDatabaseException

from Bio.Align import PairwiseAligner

class BIGSdbMLSTProfiler(AbstractAsyncContextManager):

    @abstractmethod
    def fetch_mlst_allele_variants(self, sequence_strings: Iterable[str]) -> AsyncGenerator[Allele, Any]:
        pass

    @abstractmethod
    async def fetch_mlst_st(self, alleles: Union[AsyncIterable[Allele], Iterable[Allele]]) -> MLSTProfile:
        pass

    @abstractmethod
    async def profile_string(self, sequence_strings: Iterable[str]) -> MLSTProfile:
        pass

    @abstractmethod
    def profile_multiple_strings(self, named_string_groups: AsyncIterable[Iterable[NamedString]], stop_on_fail: bool = False) -> AsyncGenerator[NamedMLSTProfile, Any]:
        pass

    @abstractmethod
    async def close(self):
        pass

class OnlineBIGSdbMLSTProfiler(BIGSdbMLSTProfiler):

    def __init__(self, database_api: str, database_name: str, schema_id: int):
        self._database_name = database_name
        self._schema_id = schema_id
        self._base_url = f"{database_api}/db/{self._database_name}/schemes/{self._schema_id}/"
        self._http_client = ClientSession(self._base_url, timeout=ClientTimeout(10000))

    async def __aenter__(self):
        return self

    async def fetch_mlst_allele_variants(self, sequence_strings: Iterable[str]) -> AsyncGenerator[Allele, Any]:
        # See https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes
        uri_path = "sequence"

        for sequence_string in sequence_strings:
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
                        partial_match_profile = PartialAllelicMatchProfile(
                            percent_identity=float(partial_match["identity"]),
                            mismatches=int(partial_match["mismatches"]),
                            gaps=int(partial_match["gaps"])
                        )
                        yield Allele(
                            allele_locus=allele_loci,
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
                allele_request_dict[allele.allele_locus].append({"allele": str(allele.allele_variant)})
        else:
            for allele in alleles:
                allele_request_dict[allele.allele_locus].append({"allele": str(allele.allele_variant)})
        request_json = {
            "designations": allele_request_dict
        }
        async with self._http_client.post(uri_path, json=request_json) as response:
            response_json: dict = await response.json()
            allele_map: dict[str, Allele] = {}
            response_json.setdefault("fields", dict())
            schema_fields_returned: dict[str, str] = response_json["fields"]
            schema_fields_returned.setdefault("ST", "unknown")
            schema_fields_returned.setdefault("clonal_complex", "unknown")
            schema_exact_matches: dict = response_json["exact_matches"]
            for exact_match_locus, exact_match_alleles in schema_exact_matches.items():
                if len(exact_match_alleles) > 1:
                    raise ValueError(f"Unexpected number of alleles returned for exact match (Expected 1, retrieved {len(exact_match_alleles)})")
                allele_map[exact_match_locus] = Allele(exact_match_locus, exact_match_alleles[0]["allele_id"], None)
            if len(allele_map) == 0:
                raise ValueError("Passed in no alleles.")
            return MLSTProfile(dict(allele_map), schema_fields_returned["ST"], schema_fields_returned["clonal_complex"])

    async def profile_string(self, sequence_strings: Iterable[str]) -> MLSTProfile:
        alleles = self.fetch_mlst_allele_variants(sequence_strings)
        return await self.fetch_mlst_st(alleles)

    async def profile_multiple_strings(self, named_string_groups: AsyncIterable[Iterable[NamedString]], stop_on_fail: bool = False) -> AsyncGenerator[NamedMLSTProfile, Any]:
        async for named_strings in named_string_groups:
            for named_string in named_strings:
                try:
                    yield NamedMLSTProfile(named_string.name, (await self.profile_string([named_string.sequence])))
                except NoBIGSdbMatchesException as e:
                    if stop_on_fail:
                        raise e
                    yield NamedMLSTProfile(named_string.name, None)

    async def close(self):
        await self._http_client.close()

    async def __aexit__(self, exc_type, exc_value, traceback):
        await self.close()

class LazyPersistentCachedBIGSdbMLSTProfiler(BIGSdbMLSTProfiler):
    def __init__(self, database_api: str, database_name: str, schema_id: int, cache_path: str):
        self._database_api = database_api
        self._database_name = database_name
        self._schema_id = schema_id
        self._base_url = f"{database_api}/db/{self._database_name}/schemes/{self._schema_id}/"
        self._http_client = ClientSession(self._base_url, timeout=ClientTimeout(10000))
        self._cache_path = cache_path
        self._loci: list[str] = []
        self._profiles = {}

    async def load_scheme_locis(self):
        self._loci.clear()
        async with self._http_client.get("") as schema_response:
            schema_json = await schema_response.json()
            for locus in schema_json["loci"]:
                locus_name = path.basename(locus)
                self._loci.append(locus_name)
        self._loci.sort()
    
    async def load_scheme_profiles(self):
        self._profiles.clear()
        with open(self.get_scheme_profile_path()) as profile_cache_handle:
            reader = csv.DictReader(profile_cache_handle, delimiter="\t")
            for line in reader:
                alleles = []
                for locus in self._loci:
                    alleles.append(line[locus])
                self._profiles[tuple(alleles)] = (line["ST"], line["clonal_complex"])
            
    def get_locus_cache_path(self, locus) -> str:
        return path.join(self._cache_path, locus + "." + "fasta")

    def get_scheme_profile_path(self):
        return path.join(self._cache_path, "profiles.csv")

    async def download_alleles_cache_data(self):
        for locus in self._loci:
            with open(self.get_locus_cache_path(locus), "wb") as fasta_handle:
                async with self._http_client.get(f"/db/{self._database_name}/loci/{locus}/alleles_fasta") as fasta_response:
                    async for chunk, eof in fasta_response.content.iter_chunks(): # TODO maybe allow chunking to be configurable
                        fasta_handle.write(chunk)

    async def download_scheme_profiles(self):
        with open(self.get_scheme_profile_path(), "wb") as profile_cache_handle:
            async with self._http_client.get("profiles_csv") as profiles_response:
                async for chunk, eof in profiles_response.content.iter_chunks():
                    profile_cache_handle.write(chunk)
    
    async def fetch_mlst_allele_variants(self, sequence_strings: Iterable[str]) -> AsyncGenerator[Allele, Any]:
        aligner = PairwiseAligner("blastn")
        aligner.mode = "local"
        for sequence_string in sequence_strings:
            for locus in self._loci:
                async for fasta_seq in read_fasta(self.get_locus_cache_path(locus)):
                    allele_variant = fasta_seq.name
                    alignment_results = aligner.align(sequence_string, fasta_seq.sequence)
                    top_alignment = sorted(alignment_results)[0]
                    top_alignment_stats = top_alignment.counts()
                    top_alignment_gaps = top_alignment_stats.gaps
                    top_alignment_identities = top_alignment_stats.identities
                    top_alignment_mismatches = top_alignment_stats.mismatches
                    if top_alignment_gaps == 0 and top_alignment_mismatches == 0:
                        yield Allele(locus, allele_variant, None)
                    else:
                        yield Allele(
                            locus,
                            allele_variant,
                            PartialAllelicMatchProfile(
                                percent_identity=top_alignment_identities/top_alignment.length,
                                mismatches=top_alignment_mismatches,
                                gaps=top_alignment_gaps
                            )
                        )

    async def fetch_mlst_st(self, alleles):
        allele_variants: dict[str, Allele] = {}
        if isinstance(alleles, AsyncIterable):
            async for allele in alleles:
                allele_variants[allele.allele_locus] = allele
        else:
            for allele in alleles:
                allele_variants[allele.allele_locus] = allele
        ordered_profile = []
        for locus in self._loci:
               ordered_profile.append(allele_variants[locus].allele_variant)

        st, clonal_complex = self._profiles[tuple(ordered_profile)]
        return MLSTProfile(allele_variants, st, clonal_complex)

    async def profile_string(self, sequence_strings: Iterable[str]) -> MLSTProfile:
        alleles = self.fetch_mlst_allele_variants(sequence_strings)
        return await self.fetch_mlst_st(alleles)

    async def profile_multiple_strings(self, named_string_groups: AsyncIterable[Iterable[NamedString]], stop_on_fail: bool = False) -> AsyncGenerator[NamedMLSTProfile, Any]:
        async for named_strings in named_string_groups:
            for named_string in named_strings:
                try:
                    yield NamedMLSTProfile(named_string.name, await self.profile_string([named_string.sequence]))
                except NoBIGSdbMatchesException as e:
                    if stop_on_fail:
                        raise e
                    yield NamedMLSTProfile(named_string.name, None)

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

    async def build_profiler_from_seqdefdb(self, dbseqdef_name: str, schema_id: int) -> OnlineBIGSdbMLSTProfiler:
        return OnlineBIGSdbMLSTProfiler(await self.get_bigsdb_api_from_seqdefdb(dbseqdef_name), dbseqdef_name, schema_id)

    async def close(self):
        await self._http_client.close()

    async def __aexit__(self, exc_type, exc_value, traceback):
        await self.close()
    
