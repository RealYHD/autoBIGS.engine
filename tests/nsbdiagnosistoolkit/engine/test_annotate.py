from automlst.engine.annotate import annotate_from_genbank, fetch_ncbi_genbank
from Bio import SeqIO

from automlst.engine.data.genomics import AnnotatedString

async def test_annotate_from_genbank_for_adk_annotation():
    sequence = str(SeqIO.read("tests/resources/tohama_I_bpertussis.fasta", "fasta").seq)
    annotated_sequence = await annotate_from_genbank("CP011448.1", "bpertussis_tohamaI", sequence, max_annotation_length=750, gene_targets=set(["adk"]))
    assert isinstance(annotated_sequence, AnnotatedString)
    assert len(annotated_sequence.annotations) >= 1
    assert annotated_sequence.annotations[0].type == "gene"
    assert "adk" in annotated_sequence.annotations[0].feature_properties["gene"]
