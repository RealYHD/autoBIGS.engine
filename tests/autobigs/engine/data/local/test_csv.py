from autobigs.engine.data.local.csv import dict_loci_alleles_variants_from_loci
from autobigs.engine.data.structures.mlst import Allele


def test_dict_loci_alleles_variants_from_loci_single_loci_not_list():
    alleles_map = {
        "adk": [Allele("adk", "1", None)]
    }
    results = dict_loci_alleles_variants_from_loci(alleles_map)
    for loci, variant in results.items():
        assert isinstance(variant, str)
        assert variant == "1"

def test_dict_loci_alleles_variants_from_loci_multi_loci_is_list():
    alleles_map = {
        "adk": [Allele("adk", "1", None), Allele("adk", "2", None)]
    }
    results = dict_loci_alleles_variants_from_loci(alleles_map)
    for loci, variant in results.items():
        assert isinstance(variant, list)
        assert len(variant) == 2