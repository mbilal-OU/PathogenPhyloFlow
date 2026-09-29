if FUNCTIONAL_ENABLED:

    rule abricate_resfinder:
        input:
            "results/input/assemblies/{sample}.fna"
        output:
            tsv="results/functional/resfinder/{sample}.tsv",
        log:
            "logs/functional/resfinder/{sample}.log"
        conda:
            "../envs/functional.yaml"
        shell:
            r"""
            mkdir -p results/functional/resfinder $(dirname {log})
            abricate --db resfinder {input} > {output.tsv} 2> {log}
            """


    rule abricate_vfdb:
        input:
            "results/input/assemblies/{sample}.fna"
        output:
            "results/functional/vfdb/{sample}.tsv"
        log:
            "logs/functional/vfdb/{sample}.log"
        conda:
            "../envs/functional.yaml"
        shell:
            "mkdir -p results/functional/vfdb $(dirname {log}) && abricate --db vfdb {input} > {output} 2> {log}"


    rule summarize_resfinder:
        input:
            expand("results/functional/resfinder/{sample}.tsv", sample=SAMPLES)
        output:
            "results/functional/resfinder_summary.tsv"
        log:
            "logs/functional/summarize_resfinder.log"
        conda:
            "../envs/functional.yaml"
        shell:
            "mkdir -p $(dirname {log}) && abricate --summary {input} > {output} 2> {log}"


    rule summarize_vfdb:
        input:
            expand("results/functional/vfdb/{sample}.tsv", sample=SAMPLES)
        output:
            "results/functional/vfdb_summary.tsv"
        log:
            "logs/functional/summarize_vfdb.log"
        conda:
            "../envs/functional.yaml"
        shell:
            "mkdir -p $(dirname {log}) && abricate --summary {input} > {output} 2> {log}"
