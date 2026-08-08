from .config import ProcessingConfig
from .rules_engine import evaluate_gene_rules


def process_paraphase_json(data: dict, config: ProcessingConfig) -> dict:
    """
    Process a single sample JSON structure, applying handlers and optionally filtering genes.
    """
    skip_keys = config.skip_keys
    genes_to_keep = (
        {g.lower() for g in config.genes_list} if config.genes_list else None
    )

    handlers = {
        "region_depth": handle_region_depth,
        "final_haplotypes": handle_final_haplotypes,
        "phase_region": handle_phase_region,
        "smn_del78_haplotypes": handle_final_haplotypes,
        "smn2_del78_haplotypes": handle_final_haplotypes,
        "smn1_haplotypes": handle_final_haplotypes,
        "smn2_haplotypes": handle_final_haplotypes,
        "fusions_called": handle_fusions_called,
        "flanking_summary": handle_dict_to_list,
    }

    if genes_to_keep:
        # Keep only selected genes (case-insensitive)
        data = {
            gene: info for gene, info in data.items() if gene.lower() in genes_to_keep
        }

    out = {}
    for gene, info in data.items():
        info = normalize_gene_info(info)
        processed = process_gene_info(info, handlers, skip_keys)

        # Optional, per-gene classification rules
        if config.rules:
            status, matches = evaluate_gene_rules(gene, processed, config.rules)
            if status is not None:
                processed["status"] = status
                # Keep a lightweight trace in json
                if matches:
                    processed["status_matches"] = [
                        {
                            "status": m.status,
                            "rule_index": m.rule_index,
                            "reason": m.rule.get("reason"),
                            "rule": m.rule,
                        }
                        for m in matches
                    ]

        out[gene] = processed

    return out


def normalize_gene_info(gene_info: dict) -> dict:
    """Flatten Paraphase 4.x region fields into the 3.x-compatible shape."""
    region_specific_info = gene_info.get("region_specific_info")
    if region_specific_info is None:
        return gene_info
    if not isinstance(region_specific_info, dict):
        raise ValueError("region_specific_info must be a mapping")

    flat_info = {
        key: value for key, value in gene_info.items() if key != "region_specific_info"
    }
    conflicting_keys = sorted(
        key
        for key in region_specific_info.keys() & flat_info.keys()
        if not _json_values_equal(region_specific_info[key], flat_info[key])
    )
    if conflicting_keys:
        raise ValueError(
            "Conflicting gene information in top-level and region_specific_info: "
            + ", ".join(conflicting_keys)
        )

    normalized = {}
    for key, value in gene_info.items():
        if key == "region_specific_info":
            normalized.update(region_specific_info)
        else:
            normalized[key] = value
    return normalized


def _json_values_equal(left, right) -> bool:
    """Compare JSON values without treating booleans as numbers."""
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(
            _json_values_equal(left[key], right[key]) for key in left
        )
    if isinstance(left, list):
        return len(left) == len(right) and all(
            _json_values_equal(left_item, right_item)
            for left_item, right_item in zip(left, right)
        )
    return left == right


def process_gene_info(gene_info, handlers, skip_keys):
    """
    Apply per-key handlers and drop skipped/None values under a gene.
    """
    processed = {}
    for gene_metric, gene_metric_value in gene_info.items():
        if gene_metric in skip_keys or gene_metric_value is None:
            continue

        if gene_metric in handlers:
            gene_metric_value = handlers[gene_metric](gene_metric_value)

        # TODO: Stringify values here instead, for both JSON and TSV output?
        processed[gene_metric] = gene_metric_value

    return processed


def handle_region_depth(value):
    """
    Replace dict, e.g. { "median": 38.0, "percentile80": 45.2 } with median value.
    """
    return value["median"]


def handle_final_haplotypes(value):
    """
    Replace a haplotypes dict , e.g.
    { "2111122111111111112111112111111111112111111111111111111111211": "smn1_smn1hap1" }
    with a list of haplotype values, e.g. ["smn1_smn1hap1"].
    """
    return list(value.values())


def handle_dict_to_list(content):
    """
    Convert dict to a list. E.g.
    "flanking_summary": {
        "f8_int22h2hap1": "region2-region2",
        "f8_int22h3hap1": "region3-region3",
        "f8_int22h1hap1": "region1-region1"
    }
    to [f8_int22h2hap1:region2-region2, f8_int22h3hap1:region3-region3, f8_int22h1hap1:region1-region1]
    """
    return [f"{key}:{value}" for key, value in content.items()]


def handle_fusions_called(content):
    """
    Remove sequence from fusions_called dict; e.g.:

        "CFH_hap2": {
            "type": "deletion",
            "sequence": "111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111122222222222222221112222122222122222222222212221112222222222222222122222222212222222222222222222222222222221222222222222222222222222222222222222222222221222222222222222222222222",
            "breakpoint": [
                [
                    19675557,
                    19760029
                ],
                [
                    19684223,
                    16844710
                ]
            ]
        }
    """
    return {
        haplotype: {
            fusion_information: fusion_value
            for fusion_information, fusion_value in haplotype_dict.items()
            if fusion_information != "sequence"
        }
        for haplotype, haplotype_dict in content.items()
    }


def handle_phase_region(content: str | None) -> list[str]:
    """
    Strip the leading prefix from each phase region, and return a list.

    Example:
        "38:chr1:196740000-196772000,38:chr1:196786000-196830000" ->
        ["chr1:196740000-196772000", "chr1:196786000-196830000"]
    """
    if not content:
        return []

    return [region.split(":", 1)[1] for region in content.split(",") if ":" in region]
