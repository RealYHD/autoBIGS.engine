import asyncio
from numbers import Number
from os import path
from typing import AsyncGenerator, Collection, Sequence, Union
from automlst.engine.data.genomics import NamedString, SangerTraceData
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Align


def _biopython_read_abif_sequence(seq_path: str) -> SeqRecord:
    with open(seq_path, "rb") as seq_handle:
        return SeqIO.read(seq_handle, "abi")


async def read_abif(seq_path: str) -> SangerTraceData:
    ext = path.splitext(seq_path)[1]
    if ext.lower() != ".ab1" and ext.lower() != "abi":
        raise ValueError(
            'seq_path must have file extension of "ab1", or "abi".')
    biopython_seq = await asyncio.to_thread(_biopython_read_abif_sequence, seq_path)
    biopython_annotations = biopython_seq.annotations

    # Lot of type ignoring since Biopython did not define their typing.
    biopython_abif_raw = biopython_annotations["abif_raw"] # type: ignore
    trace_data = SangerTraceData(
        path.basename(seq_path),
        biopython_seq.seq,
        biopython_abif_raw.get("APFN2"), # type: ignore
        biopython_abif_raw.get("APrN1"), # type: ignore
        biopython_abif_raw.get("APrV1"), # type: ignore
        biopython_abif_raw.get("APrX1"), # type: ignore
        biopython_abif_raw.get("APXV1"), # type: ignore
        biopython_abif_raw.get("CMNT1"), # type: ignore
        biopython_abif_raw.get("CpEP1"), # type: ignore
        biopython_abif_raw.get("CTID1"), # type: ignore
        biopython_abif_raw.get("CTNM1"), # type: ignore
        biopython_abif_raw.get("CTTL1"), # type: ignore
        biopython_abif_raw.get("DATA1"), # type: ignore
        biopython_abif_raw.get("DATA2"), # type: ignore
        biopython_abif_raw.get("DATA3"), # type: ignore
        biopython_abif_raw.get("DATA4"), # type: ignore
        biopython_abif_raw.get("DATA5"), # type: ignore
        biopython_abif_raw.get("DATA6"), # type: ignore
        biopython_abif_raw.get("DATA7"), # type: ignore
        biopython_abif_raw.get("DATA8"), # type: ignore
        biopython_abif_raw.get("DSam1"), # type: ignore
        biopython_abif_raw.get("DyeN1"), # type: ignore
        biopython_abif_raw.get("DyeN2"), # type: ignore
        biopython_abif_raw.get("DyeN3"), # type: ignore
        biopython_abif_raw.get("DyeN4"), # type: ignore
        biopython_abif_raw.get("DyeW1"), # type: ignore
        biopython_abif_raw.get("DyeW2"), # type: ignore
        biopython_abif_raw.get("DyeW3"), # type: ignore
        biopython_abif_raw.get("DyeW4"), # type: ignore
        biopython_abif_raw.get("DySN1"), # type: ignore
        biopython_abif_raw.get("EPVt1"), # type: ignore
        biopython_abif_raw.get("EVNT1"), # type: ignore
        biopython_abif_raw.get("EVNT2"), # type: ignore
        biopython_abif_raw.get("EVNT3"), # type: ignore
        biopython_abif_raw.get("EVNT4"), # type: ignore
        biopython_abif_raw.get("FWO_1"), # type: ignore
        biopython_abif_raw.get("GTyp1"), # type: ignore
        biopython_abif_raw.get("InSc1"), # type: ignore
        biopython_abif_raw.get("InVt1"), # type: ignore
        biopython_abif_raw.get("LANE1"), # type: ignore
        biopython_abif_raw.get("LIMS1"), # type: ignore
        biopython_abif_raw.get("LNTD1"), # type: ignore
        biopython_abif_raw.get("LsrP1"), # type: ignore
        biopython_abif_raw.get("MCHN1"), # type: ignore
        biopython_abif_raw.get("MODF1"), # type: ignore
        biopython_abif_raw.get("MODL1"), # type: ignore
        biopython_abif_raw.get("NAVG1"), # type: ignore
        biopython_abif_raw.get("NLNE1"), # type: ignore
        biopython_abif_raw.get("OfSc1"), # type: ignore
        biopython_abif_raw.get("PDMF1"), # type: ignore
        biopython_abif_raw.get("PXLB1"), # type: ignore
        biopython_abif_raw.get("RGCm1"), # type: ignore
        biopython_abif_raw.get("RGNm1"), # type: ignore
        biopython_abif_raw.get("RMdV1"), # type: ignore
        biopython_abif_raw.get("RMdX1"), # type: ignore
        biopython_abif_raw.get("RMXV1"), # type: ignore
        biopython_abif_raw.get("RPrN1"), # type: ignore
        biopython_abif_raw.get("RPrV1"), # type: ignore
        biopython_abif_raw.get("RUND1"), # type: ignore
        biopython_abif_raw.get("RUND2"), # type: ignore
        biopython_abif_raw.get("RUND3"), # type: ignore
        biopython_abif_raw.get("RUND4"), # type: ignore
        biopython_abif_raw.get("RunN1"), # type: ignore
        biopython_abif_raw.get("RUNT1"), # type: ignore
        biopython_abif_raw.get("RUNT2"), # type: ignore
        biopython_abif_raw.get("RUNT3"), # type: ignore
        biopython_abif_raw.get("RUNT4"), # type: ignore
        biopython_abif_raw.get("Satd"), # type: ignore
        biopython_abif_raw.get("Scal1"), # type: ignore
        biopython_abif_raw.get("SCAN1"), # type: ignore
        biopython_abif_raw.get("SMED1"), # type: ignore
        biopython_abif_raw.get("SMLt"), # type: ignore
        biopython_abif_raw.get("SMPL1"), # type: ignore
        biopython_abif_raw.get("SVER1"), # type: ignore
        biopython_abif_raw.get("SVER3"), # type: ignore
        biopython_abif_raw.get("Tmpr1"), # type: ignore
        biopython_abif_raw.get("TUBE"), # type: ignore
        biopython_abif_raw.get("User") # type: ignore
    )
    return trace_data

def _biopython_local_pairwise_alignment(reference: NamedString, query: NamedString) -> tuple[NamedString, NamedString]:
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "local"
    alignment_result = sorted(aligner.align(reference.sequence, query.sequence))[0] # take the best alignment
    return NamedString(alignment_result.sequences[0].id, alignment_result.sequences[0].seq), NamedString(alignment_result.sequences[1].id, alignment_result.sequences[1].seq)

async def reference_consensus_assembly(reference: NamedString, sanger_traces: Collection[SangerTraceData]) -> AsyncGenerator[NamedString, Any]:
    for sanger_trace in sanger_traces:
        yield (await asyncio.to_thread(_biopython_local_pairwise_alignment, reference, sanger_trace))[1]