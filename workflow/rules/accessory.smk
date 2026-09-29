import os

if ACCESSORY_ENABLED:

    rule prokka_sample:
        input:
            "results/input/assemblies/{sample}.fna"
        output:
            gff="results/accessory/prokka/{sample}/{sample}.gff"
        params:
            # GNU parallel --pipe workers exec the shell named by the SHELL
            # environment variable, which is not guaranteed to be exported in
            # non-interactive job environments (snakemake, cron, containers);
            # without it parallel dies with "exec failed: No such file or
            # directory". PARALLEL_SHELL is honored by parallel instead, so we
            # pass an explicit fallback here and still respect a user-provided
            # PARALLEL_SHELL when one is set.
            parallel_shell=lambda wildcards: os.environ.get(
                "PARALLEL_SHELL", "/bin/bash"
            ),
        threads:
            config["resources"]["prokka_threads"]
        resources:
            mem_mb=config["resources"].get("prokka_mem_mb", 4000),
        log:
            "logs/accessory/prokka/{sample}.log"
        conda:
            "../envs/prokka.yaml"
        shell:
            r"""
            mkdir -p results/accessory/prokka/{wildcards.sample} $(dirname {log})
            PARALLEL_SHELL={params.parallel_shell} \
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
            # Panaroo's --alignment mode. "none" skips gene alignment entirely
            # (much faster); the core/pan alignments are not consumed by any
            # downstream rule, which only needs the presence/absence matrices.
            # Default "core" preserves the original behavior.
            alignment=config["accessory"].get("panaroo_alignment", "core"),
            gffs=" ".join(
                f"results/accessory/prokka/{sample}/{sample}.gff" for sample in SAMPLES
            ),
        threads:
            config["resources"]["panaroo_threads"]
        resources:
            mem_mb=config["resources"].get("panaroo_mem_mb", 16000),
        log:
            "logs/accessory/panaroo.log"
        conda:
            "../envs/panaroo.yaml"
        shell:
            r"""
            mkdir -p results/accessory/panaroo $(dirname {log})
            ALIGN=""
            if [ "{params.alignment}" != "none" ]; then
                ALIGN="--alignment {params.alignment}"
            fi
            panaroo -i {params.gffs} -o results/accessory/panaroo \
              --clean-mode {params.clean} --core_threshold {params.core} \
              $ALIGN --threads {threads} --remove-invalid-genes \
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
