rule select_reference:
    input:
        assemblies=expand("results/input/assemblies/{sample}.fna", sample=SAMPLES)
    output:
        reference="results/reference/selected_reference.fna",
        table="results/reference/reference_candidates.tsv",
        summary="results/reference/reference_selection.json",
    params:
        samples=SAMPLES,
        settings=config["reference"],
    log:
        "logs/reference/select_reference.log"
    conda:
        "../envs/reference.yaml"
    script:
        "../scripts/reference_select.py"
