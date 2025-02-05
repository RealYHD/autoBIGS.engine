import asyncio
from concurrent.futures import Future, ThreadPoolExecutor
from contextlib import AbstractContextManager
from typing import Any, Set, Union
from Bio.Align import PairwiseAligner
from queue import Queue

from autobigs.engine.structures.alignment import AlignmentStats, PairwiseAlignment


class AsyncPairwiseAlignmentEngine(AbstractContextManager):
    def __enter__(self):
        self._thread_pool = ThreadPoolExecutor(self._max_threads, thread_name_prefix="async-pairwise-alignment")
        return self

    def __init__(self, aligner: PairwiseAligner, max_threads: int = 4):
        self._max_threads = max_threads
        self._aligner = aligner
        self._work_left: Set[Future] = set()
        self._work_complete: Queue[Future] = Queue()

    def align(self, reference: str, query: str, **associated_data):
        work = self._thread_pool.submit(
            self.work, reference, query, **associated_data)
        work.add_done_callback(self._on_complete)
        self._work_left.add(work)
        
    def _on_complete(self, future: Future):
        self._work_left.remove(future)
        self._work_complete.put(future)

    def work(self, reference, query, **associated_data):
        alignments = self._aligner.align(reference, query)
        top_alignment = alignments[0]
        top_alignment_stats = top_alignment.counts()
        top_alignment_gaps = top_alignment_stats.gaps
        top_alignment_identities = top_alignment_stats.identities
        top_alignment_mismatches = top_alignment_stats.mismatches
        top_alignment_score = top_alignment.score # type: ignore
        return PairwiseAlignment(
            top_alignment.sequences[0],
            top_alignment.sequences[1],
            tuple(top_alignment.indices[0]),
            tuple(top_alignment.indices[1]),
            AlignmentStats(
                percent_identity=top_alignment_identities/top_alignment.length,
                mismatches=top_alignment_mismatches,
                gaps=top_alignment_gaps,
                score=top_alignment_score
            )), associated_data

    async def next_completed(self) -> Union[tuple[PairwiseAlignment, dict[str, Any]], None]:
        if self._work_complete.empty() and len(self._work_left):
            return None
        future_now = await asyncio.wrap_future(self._work_complete.get())
        return future_now

    def __exit__(self, exc_type, exc_value, traceback):
        self.shutdown()

    def __aiter__(self):
        return self
    
    async def __anext__(self):
        result = await self.next_completed()
        if result is None:
            raise StopAsyncIteration
        return result

    def shutdown(self):
        self._thread_pool.shutdown(wait=True, cancel_futures=True)
