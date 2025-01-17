from automlst.engine.data.local.csv import dict_loci_alleles_variants_from_loci
from automlst.engine.data.structures.mlst import Allele


def test_dict_loci_alleles_variants_from_loci_single_loci_not_list():
    alleles_map = {
        "adk": [Allele("adk", "1", None)]
    }
    results = dict_loci_alleles_variants_from_loci(alleles_map)
    for loci, variant in results.items():
        assert isinstance(variant, str)
        assert variant == "1"