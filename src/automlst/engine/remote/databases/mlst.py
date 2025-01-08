from abc import abstractmethod
from contextlib import AbstractAsyncContextManager
from typing import AsyncGenerator, AsyncIterable, Generator, Iterable, Mapping, Union

from aiohttp import ClientSession

from automlst.engine.data.MLST import Allele, MLSTProfile

MLST_DATABASES = [
    "https://bigsdb.pasteur.fr/api/db",
    "https://rest.pubmlst.org/db"
]

class MLSTProfiler(AbstractAsyncContextManager):
    @abstractmethod
    def fetch_mlst_allele_variants(self, schema_id: int, sequence_string: str) -> AsyncGenerator[Allele]:
        pass
    
    @abstractmethod
    async def fetch_mlst_st(self, schema_id: int, alleles: AsyncIterable[Allele]) -> MLSTProfile:
        pass

    @abstractmethod
    async def profile_string(self, schema_id: int, string: str) -> MLSTProfile:
        pass

    @abstractmethod
    async def close(self):
        pass

    @abstractmethod
    async def get_scheme_ids(self) -> Mapping[str, int]:
        pass