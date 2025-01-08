import asyncio
from Bio import Entrez
from Bio import SeqIO

# TODO Change this out for a more professional approach
Entrez.email = "yunyangdeng@outlook.com"

from automlst.engine.data.genomics import AnnotatedString, StringAnnotation


async def fetch_ncbi_genbank(genbank_id: str) -> AnnotatedString:
    with (await asyncio.to_thread(Entrez.efetch, db="nucleotide", id=genbank_id, rettype="gb", retmode="text")) as fetch_stream:
        record = SeqIO.read(fetch_stream, "genbank")
        sequence_features = list()
        for feature in record.features:
            start = int(feature.location.start)
            end = int(feature.location.end)
            qualifiers = feature.qualifiers
            for qualifier_key in qualifiers:
                qualifiers[qualifier_key] = set(qualifiers[qualifier_key])
            sequence_features.append(StringAnnotation(
                type=feature.type,
                start=start,
                end=end+1,  # Position is exclusive
                feature_properties=qualifiers
            ))
        return AnnotatedString(name=genbank_id, sequence=str(record.seq), annotations=sequence_features)