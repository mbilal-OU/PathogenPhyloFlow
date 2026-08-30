if ACCESSORY_ENABLED:

    rule prokka_sample:
        input:
            "results/input/assemblies/{sample}.fna"
        output:
            gff="results/accessory/prokka/{sample}/{sample}.gff"
        threads:
            config["resources"]["prokka_threads"]
        log:
            "logs/accessory/prokka/{sample}.log"
        conda:
            "../envs/prokka.yaml"
        shell:
            r"""
            mkdir -p results/accessory/prokka/{wildcards.sample} $(dirname {log})
            prokka --outdir results/accessory/prokka/{wildcards.sample} \
              --prefix {wildcards.sample} --locustag {wildcards.sample} \
              --cpus {threads} --force {input} > {log} 2>&1
            """


    rule panaroo:
        input:
            gffs=expand("results/accessory/prokka/{sample}/{sample}.gff", sample=SAMPLES)
        output:
            rtab="results/accessory/panaroo/gene_presence_absence.Rtab",
            csv="results/accessory/panaroo/gene_presence_absence.csv",
        params:
            clean=config["accessory"]["clean_mode"],
            core=config["accessory"]["core_threshold"],
            gffs=" ".join(
                f"results/accessory/prokka/{sample}/{sample}.gff" for sample in SAMPLES
            ),
        threads:
            config["resources"]["panaroo_threads"]
        log:
            "logs/accessory/panaroo.log"
        conda:
            "../envs/panaroo.yaml"
        shell:
            r"""
            mkdir -p results/accessory/panaroo $(dirname {log})
            panaroo -i {params.gffs} -o results/accessory/panaroo \
              --clean-mode {params.clean} --core_threshold {params.core} \
              --alignment core --threads {threads} --remove-invalid-genes \
              > {log} 2>&1
            """


    rule snp_accessory_discordance:
        input:
            snp="results/phylogeny/final.snp_dist.tsv",
            accessory="results/accessory/panaroo/gene_presence_absence.Rtab",
        output:
            table="results/integration/snp_accessory_discordance.tsv",
            summary="results/integration/snp_accessory_summary.json",
        params:
            low_snp=config["integration"]["low_snp_threshold"],
            high_accessory=config["integration"]["high_accessory_distance"],
        conda:
            "../envs/report.yaml"
        script:
            "../scripts/discordance.py"
