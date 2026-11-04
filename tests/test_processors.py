"""Basic tests for paraphrase."""


def test_process_paraphase_json():
    from paraphrase.processors import process_paraphase_json
    from paraphrase.config import ProcessingConfig

    input_data = {
        "smn1": {
            "region_depth": {"median": 44.0, "percentile80": 48.0},
            "final_haplotypes": {
                "x22222122222221111111111": "smn1_smn2hap1",
                "x22222222221112222222222": "smn1_smn2hap2",
            },
        },
        "F8": {
            "flanking_summary": {
                "f8_int22h2hap1": "region2-region2",
                "f8_int22h3hap1": "-region3",
                "f8_int22h2hap2": "-region2",
                "f8_int22h1hap1": "region1-region1",
                "f8_int22h1hap2": "region1-region1",
            }
        },
    }

    config = ProcessingConfig(
        skip_keys=set(),
        genes_list=["smn1", "f8"],
        rules=None,
    )

    expected_output = {
        "smn1": {
            "region_depth": 44.0,
            "final_haplotypes": ["smn1_smn2hap1", "smn1_smn2hap2"],
        },
        "F8": {
            "flanking_summary": [
                "f8_int22h2hap1:region2-region2",
                "f8_int22h3hap1:-region3",
                "f8_int22h2hap2:-region2",
                "f8_int22h1hap1:region1-region1",
                "f8_int22h1hap2:region1-region1",
            ]
        },
    }
    assert process_paraphase_json(input_data, config) == expected_output


def test_process_gene_info():
    from paraphrase.processors import process_gene_info
    from paraphrase.processors import (
        handle_region_depth,
        handle_final_haplotypes,
    )

    handlers = {
        "region_depth": handle_region_depth,
        "final_haplotypes": handle_final_haplotypes,
    }
    gene_info = {
        "region_depth": {"median": 44.0, "percentile80": 48.0},
        "final_haplotypes": {
            "x22222122222221111111111": "smn1_smn2hap1",
            "x22222222221112222222222": "smn1_smn2hap2",
        },
    }
    skip_keys = {"final_haplotypes"}
    expected_output = {"region_depth": 44.0}
    assert process_gene_info(gene_info, handlers, skip_keys) == expected_output

    # The remaining tests in this module exercised the old `normal_values`-based
    # per-key flagging logic, which has been superseded by the rules engine and
    # removed from `process_gene_info`. They are intentionally omitted here.


def test_handle_region_depth():
    from paraphrase.processors import handle_region_depth

    input_data = {"median": 44.0, "percentile80": 48.0}

    expected_output = 44.0
    assert handle_region_depth(input_data) == expected_output


def test_process_gene_info_drops_empty_region_depth():
    """Paraphase 4.x reports region_depth as {} for some regions (e.g. CFHclust)."""
    from paraphrase.processors import process_gene_info, handle_region_depth

    handlers = {"region_depth": handle_region_depth}
    gene_info = {"region_depth": {}, "total_cn": 3}
    assert process_gene_info(gene_info, handlers, set()) == {"total_cn": 3}


def test_handle_final_haplotypes():
    from paraphrase.processors import handle_final_haplotypes

    input_data = {
        "x22222122222221111111111": "smn1_smn2hap1",
        "x22222222221112222222222": "smn1_smn2hap2",
    }

    expected_output = ["smn1_smn2hap1", "smn1_smn2hap2"]
    assert handle_final_haplotypes(input_data) == expected_output


def test_handle_dict_to_list():
    from paraphrase.processors import handle_dict_to_list

    input_data = {
        "f8_int22h2hap1": "region2-region2",
        "f8_int22h3hap1": "-region3",
        "f8_int22h2hap2": "-region2",
        "f8_int22h1hap1": "region1-region1",
        "f8_int22h1hap2": "region1-region1",
    }

    expected_output = [
        "f8_int22h2hap1:region2-region2",
        "f8_int22h3hap1:-region3",
        "f8_int22h2hap2:-region2",
        "f8_int22h1hap1:region1-region1",
        "f8_int22h1hap2:region1-region1",
    ]
    assert handle_dict_to_list(input_data) == expected_output


def test_handle_fusions_called():
    from paraphrase.processors import handle_fusions_called

    input_data = {
        "CFH_hap1": {
            "type": "deletion",
            "sequence": "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111122222222222222221112222122222222222222222222221112222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222",
            "breakpoint": [[196757557, 196760029], [196842234, 196844710]],
        }
    }

    expected_output = {
        "CFH_hap1": {
            "type": "deletion",
            "breakpoint": [[196757557, 196760029], [196842234, 196844710]],
        }
    }
    assert handle_fusions_called(input_data) == expected_output
