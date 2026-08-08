"""Basic tests for paraphrase."""

import pytest


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


def test_paraphase_v4_region_fields_match_v3_processing():
    from paraphrase.config import ProcessingConfig
    from paraphrase.processors import process_paraphase_json

    fusion = {
        "cfh_hap1": {
            "type": "deletion",
            "sequence": "121122",
            "breakpoint": [[1, 2], [3, 4]],
        }
    }
    v3 = {
        "cfh": {"total_cn": 4, "fusions_called": fusion},
        "ncf1": {"total_cn": 4, "gene_cn": 0, "gene_reads": 20},
        "strc": {"total_cn": 4, "gene_cn": None},
    }
    v4 = {
        "cfh": {
            "total_cn": 4,
            "region_specific_info": {"fusions_called": fusion},
        },
        "ncf1": {
            "total_cn": 4,
            "region_specific_info": {"gene_cn": 0, "gene_reads": 20},
        },
        "strc": {
            "total_cn": 4,
            "region_specific_info": {"gene_cn": None},
        },
    }
    rules = {
        "cfh": {
            "rules": [
                {
                    "status": "intermediate",
                    "when": {"fusions_called": {"not_empty": True}},
                }
            ]
        },
        "ncf1": {
            "rules": [
                {
                    "status": "pathological",
                    "when": {"gene_cn": 0},
                }
            ]
        },
    }
    config = ProcessingConfig(skip_keys={"gene_reads"}, rules=rules)

    v3_output = process_paraphase_json(v3, config)
    v4_output = process_paraphase_json(v4, config)

    assert v4_output == v3_output
    assert "sequence" not in v4_output["cfh"]["fusions_called"]["cfh_hap1"]
    assert v4_output["cfh"]["status"] == "intermediate"
    assert v4_output["ncf1"]["status"] == "pathological"
    assert "gene_reads" not in v4_output["ncf1"]
    assert "gene_cn" not in v4_output["strc"]


@pytest.mark.parametrize(
    ("flat", "nested"),
    [(2, 3), (1, True), ({"value": 1}, {"value": True}), ([1], [True])],
)
def test_process_paraphase_rejects_conflicting_flat_and_nested_values(flat, nested):
    from paraphrase.config import ProcessingConfig
    from paraphrase.processors import process_paraphase_json

    input_data = {
        "strc": {
            "gene_cn": flat,
            "region_specific_info": {"gene_cn": nested},
        }
    }

    with pytest.raises(ValueError, match="Conflicting gene information.*gene_cn"):
        process_paraphase_json(input_data, ProcessingConfig(skip_keys=set()))


def test_process_paraphase_accepts_equal_duplicate_values():
    from paraphrase.config import ProcessingConfig
    from paraphrase.processors import process_paraphase_json

    input_data = {
        "strc": {
            "gene_cn": 2,
            "region_specific_info": {"gene_cn": 2},
        }
    }

    assert process_paraphase_json(input_data, ProcessingConfig(skip_keys=set())) == {
        "strc": {"gene_cn": 2}
    }


def test_process_paraphase_expands_region_fields_in_place():
    from paraphrase.config import ProcessingConfig
    from paraphrase.processors import process_paraphase_json

    output = process_paraphase_json(
        {
            "strc": {
                "before": 1,
                "region_specific_info": {"gene_cn": 2},
                "after": 3,
            }
        },
        ProcessingConfig(skip_keys=set()),
    )

    assert list(output["strc"]) == ["before", "gene_cn", "after"]


def test_process_paraphase_rejects_non_mapping_region_specific_info():
    from paraphrase.config import ProcessingConfig
    from paraphrase.processors import process_paraphase_json

    with pytest.raises(ValueError, match="region_specific_info must be a mapping"):
        process_paraphase_json(
            {"strc": {"region_specific_info": []}},
            ProcessingConfig(skip_keys=set()),
        )


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
