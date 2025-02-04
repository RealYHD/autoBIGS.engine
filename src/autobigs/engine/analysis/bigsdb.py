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

from autobigs.engine.analysis.aligners import AsyncPairwiseAlignmentEngine
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
        self._http_client = ClientSession(self._base_url, timeout=ClientTimeout(10000))

    async def __aenter__(self):
        return self

    async def determine_mlst_allele_variants(self, query_sequence_strings: Iterable[str]) -> AsyncGenerator[Allele, Any]:
        # See https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes
        uri_path = "sequence"
        if not isinstance(query_sequence_strings, Iterable):
            raise ValueError("Invalid data type for parameter \"sequence_strings\".")

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
                            score=int(partial_match["score"])
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

class LocalBIGSdbMLSTProfiler(BIGSdbMLSTProfiler):
    async def __aenter__(self):
        if self._prepare:
            await self.update_scheme_locis()
            await asyncio.gather(
                self.download_alleles_cache_data(),
                self.download_scheme_profiles()
            )
            await self.load_scheme_profiles()
        return self
    
    def __init__(self, database_api: str, database_name: str, schema_id: int, cache_path: Union[str, None] = None, prepare: bool =True):
        self._database_api = database_api
        self._database_name = database_name
        self._schema_id = schema_id
        self._base_url = f"{self._database_api}/db/{self._database_name}/schemes/{self._schema_id}/"
        self._http_client = ClientSession(self._base_url, timeout=ClientTimeout(10000))
        if cache_path is None:
            self._cache_path = tempfile.mkdtemp("BIGSdb")
            self._cleanup_required = True
        else:
            self._cache_path = cache_path
            self._cleanup_required = False
        self._loci: list[str] = []
        self._profiles_st_map = {}
        self._prepare = prepare

    async def update_scheme_locis(self):
        self._loci.clear()
        async with self._http_client.get(f"/api/db/{self._database_name}/schemes/{self._schema_id}") as schema_response:
            schema_json = await schema_response.json()
            for locus in schema_json["loci"]:
                locus_name = path.basename(locus)
                self._loci.append(locus_name)
        self._loci.sort()
    
    async def load_scheme_profiles(self):
        self._profiles_st_map.clear()
        with open(self.get_scheme_profile_path()) as profile_cache_handle:
            reader = csv.DictReader(profile_cache_handle, delimiter="\t")
            for line in reader:
                alleles = []
                for locus in self._loci:
                    alleles.append(line[locus])
                self._profiles_st_map[tuple(alleles)] = (line["ST"], line["clonal_complex"])
            
    def get_locus_cache_path(self, locus) -> str:
        return path.join(self._cache_path, locus + "." + "fasta")

    def get_scheme_profile_path(self):
        return path.join(self._cache_path, "profiles.csv")

    async def download_alleles_cache_data(self):
        for locus in self._loci:
            with open(self.get_locus_cache_path(locus), "wb") as fasta_handle:
                async with self._http_client.get(f"/api/db/{self._database_name}/loci/{locus}/alleles_fasta") as fasta_response:
                    async for chunk, eof in fasta_response.content.iter_chunks():
                        fasta_handle.write(chunk)

    async def download_scheme_profiles(self):
        with open(self.get_scheme_profile_path(), "wb") as profile_cache_handle:
            async with self._http_client.get("profiles_csv") as profiles_response:
                async for chunk, eof in profiles_response.content.iter_chunks():
                    profile_cache_handle.write(chunk)
        await self.load_scheme_profiles()
    
    async def determine_mlst_allele_variants(self, query_sequence_strings: Iterable[str]) -> AsyncGenerator[Allele, Any]:
        aligner = PairwiseAligner("blastn")
        aligner.mode = "local"
        with AsyncPairwiseAlignmentEngine(aligner) as aligner_engine:
            for query_sequence_string in query_sequence_strings:
                for locus in self._loci:
                    async for allele_variant in read_fasta(self.get_locus_cache_path(locus)):
                        aligner_engine.align(allele_variant.sequence, query_sequence_string, variant_name=allele_variant.name, full=True)
                        break # start a bunch of full alignments for each variant to select segments
            alignment_rankings: dict[str, set[tuple[PairwiseAlignment, str]]] = defaultdict(set)
            async for alignment_result, additional_information in aligner_engine:
                result_variant_name = additional_information["variant_name"]
                result_locus, variant_id = result_variant_name.split("_")
                full_alignment = additional_information["full"]
                if full_alignment:
                    if alignment_result.alignment_stats.gaps == 0 and alignment_result.alignment_stats.mismatches == 0:
                        # I.e., 100% exactly the same
                        yield Allele(result_locus, variant_id, None)
                        continue
                    else:
                        alignment_rankings[result_locus].add((alignment_result, variant_id))
                    interest_sequence = full_alignment[alignment_result.query_indices[0]:alignment_result.query_indices[-1]]
                    async for allele_variant in read_fasta(self.get_locus_cache_path(result_locus)):
                        if result_variant_name == allele_variant.name:
                            continue # Skip if we just finished aligning this
                        aligner_engine.align(allele_variant.sequence, interest_sequence, variant_name=result_variant_name.name, full=False)
                else:
                    alignment_rankings[result_locus].add((alignment_result, variant_id))
            for final_locus, alignments in alignment_rankings.items():
                closest_alignment, closest_variant_id = sorted(alignments, key=lambda index: index[0].alignment_stats.score)[0]
                yield Allele(final_locus, closest_variant_id, closest_alignment.alignment_stats)

    async def determine_mlst_st(self, alleles):
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

        st, clonal_complex = self._profiles_st_map[tuple(ordered_profile)]
        return MLSTProfile(set(allele_variants.values()), st, clonal_complex)

    async def profile_string(self, query_sequence_strings: Iterable[str]) -> MLSTProfile:
        alleles = self.determine_mlst_allele_variants(query_sequence_strings)
        return await self.determine_mlst_st(alleles)

    async def profile_multiple_strings(self, query_named_string_groups: AsyncIterable[Iterable[NamedString]], stop_on_fail: bool = False) -> AsyncGenerator[NamedMLSTProfile, Any]:
        async for named_strings in query_named_string_groups:
            for named_string in named_strings:
                try:
                    yield NamedMLSTProfile(named_string.name, await self.profile_string([named_string.sequence]))
                except NoBIGSdbMatchesException as e:
                    if stop_on_fail:
                        raise e
                    yield NamedMLSTProfile(named_string.name, None)

    async def close(self):
        await self._http_client.close()
        if self._cleanup_required:
            shutil.rmtree(self._cache_path)

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

    async def build_profiler_from_seqdefdb(self, dbseqdef_name: str, schema_id: int) -> RemoteBIGSdbMLSTProfiler:
        return RemoteBIGSdbMLSTProfiler(await self.get_bigsdb_api_from_seqdefdb(dbseqdef_name), dbseqdef_name, schema_id)

    async def close(self):
        await self._http_client.close()

    async def __aexit__(self, exc_type, exc_value, traceback):
        await self.close()

def get_BIGSdb_MLST_profiler(local: bool, database_api: str, database_name: str, schema_id: int):
    if local:
        return LocalBIGSdbMLSTProfiler(database_api=database_api, database_name=database_name, schema_id=schema_id)
    return RemoteBIGSdbMLSTProfiler(database_api=database_api, database_name=database_name, schema_id=schema_id)