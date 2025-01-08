import asyncio
from collections.abc import Set
from typing import Any, Generator, List, Sequence
from Bio.Align import PairwiseAligner
from Bio import Entrez
from Bio import SeqIO
import numpy as np

from automlst.engine.data.genomics import AnnotatedString, StringAnnotation
from automlst.engine.remote.databases.ncbi.genbank import fetch_ncbi_genbank


async def annotate_from_genbank(genbank_id: str, query_name: str, query_string: str, max_annotation_length:int = 512, gene_targets:Set = set()):
    # TODO implement asynchronous alignment algorithm
    reference_annotations = await fetch_ncbi_genbank(genbank_id=genbank_id)
    query_annotations = list()
    aligner = PairwiseAligner("blastn")
    aligner.mode = "local"
    for annotation in reference_annotations.annotations:
        if annotation.type != "gene" or "gene" not in annotation.feature_properties:
            continue
        if len(gene_targets) > 0 and "gene" in annotation.feature_properties:
            if not annotation.feature_properties["gene"].intersection(gene_targets):
                continue
        if max_annotation_length > 0 and annotation.end - annotation.start > max_annotation_length:
            # TODO implement a failsafe
            continue
        feature_string_sequence = get_feature_coding(annotated_string=reference_annotations, string_annotation=annotation)
        alignments = aligner.align(query_string, feature_string_sequence)
        if len(alignments) < 1:
            # TODO implement a failsafe
            continue
        top_alignment = sorted(alignments)[0]
        # TODO Check if alternatives are better
        query_annotations.append(StringAnnotation(
            type=annotation.type, # same as original
            start=np.min(top_alignment.aligned[0]), # We only care about the start of first chunk
            end=np.max(top_alignment.aligned[0]), # and the end of the last chunk
            feature_properties=dict(annotation.feature_properties) # same as original
        ))
    return AnnotatedString(name=query_name, sequence=query_string, annotations=query_annotations)

def get_feature_coding(annotated_string: AnnotatedString, string_annotation: StringAnnotation) -> str:
    return annotated_string.sequence[string_annotation.start:string_annotation.end]