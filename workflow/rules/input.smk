rule resolve_input:
    output:
        assembly="results/input/assemblies/{sample}.fna",
        metadata="results/input/metadata/{sample}.json",
    params:
        row=lambda wildcards: SAMPLE_META[wildcards.sample]
    log:
        "logs/input/{sample}.log"
    conda:
        "../envs/input.yaml"
    script:
        "../scripts/resolve_input.py"
