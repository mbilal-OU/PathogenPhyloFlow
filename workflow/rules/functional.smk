if FUNCTIONAL_ENABLED:

    rule abricate_resfinder:
        input:
            "results/input/assemblies/{sample}.fna"
        output:
            "results/functional/resfinder/{sample}.tsv"
        conda:
            "../envs/functional.yaml"
        shell:
            "mkdir -p results/functional/resfinder && abricate --db resfinder {input} > {output}"


    rule abricate_vfdb:
        input:
            "results/input/assemblies/{sample}.fna"
        output:
            "results/functional/vfdb/{sample}.tsv"
        conda:
            "../envs/functional.yaml"
        shell:
            "mkdir -p results/functional/vfdb && abricate --db vfdb {input} > {output}"


    rule summarize_resfinder:
        input:
            expand("results/functional/resfinder/{sample}.tsv", sample=SAMPLES)
        output:
            "results/functional/resfinder_summary.tsv"
        conda:
            "../envs/functional.yaml"
        shell:
            "abricate --summary {input} > {output}"


    rule summarize_vfdb:
        input:
            expand("results/functional/vfdb/{sample}.tsv", sample=SAMPLES)
        output:
            "results/functional/vfdb_summary.tsv"
        conda:
            "../envs/functional.yaml"
        shell:
            "abricate --summary {input} > {output}"
