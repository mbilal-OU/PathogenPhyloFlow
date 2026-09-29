rule snippy_sample:
    input:
        assembly="results/input/assemblies/{sample}.fna",
        reference="results/reference/selected_reference.fna",
    output:
        vcf="results/variants/snippy/{sample}/snps.vcf",
    threads:
        config["resources"]["snippy_threads"]
    resources:
        # Keep the scheduler-visible memory in sync with --ram above.
        mem_mb=config["resources"].get("snippy_ram_gb", 8) * 1000,
    params:
        # Snippy defaults to 8 GB RAM for samtools sort; make it configurable
        # so low-memory machines can run it. Default 8 preserves old behavior.
        ram=config["resources"].get("snippy_ram_gb", 8),
    log:
        "logs/snippy/{sample}.log"
    conda:
        "../envs/snippy.yaml"
    shell:
        r"""
        mkdir -p results/variants/snippy/{wildcards.sample} $(dirname {log})
        snippy \
          --outdir results/variants/snippy/{wildcards.sample} \
          --ref {input.reference} \
          --ctgs {input.assembly} \
          --cpus {threads} \
          --ram {params.ram} \
          --force > {log} 2>&1
        test -s {output.vcf}
        """


rule snippy_core:
    input:
        vcfs=expand("results/variants/snippy/{sample}/snps.vcf", sample=SAMPLES),
        reference="results/reference/selected_reference.fna",
    output:
        aln="results/variants/core/core.aln",
        full="results/variants/core/core.full.aln",
        tab="results/variants/core/core.tab",
    params:
        dirs=" ".join(f"results/variants/snippy/{sample}" for sample in SAMPLES)
    log:
        "logs/snippy/core.log"
    conda:
        "../envs/snippy.yaml"
    shell:
        r"""
        mkdir -p results/variants/core $(dirname {log})
        snippy-core --ref {input.reference} --prefix results/variants/core/core {params.dirs} > {log} 2>&1
        test -s {output.aln}
        test -s {output.full}
        """


rule raw_snp_distances:
    input:
        "results/variants/core/core.aln"
    output:
        "results/variants/core/raw.snp_dist.tsv"
    conda:
        "../envs/phylogeny.yaml"
    shell:
        "snp-dists {input} > {output}"


rule raw_iqtree:
    input:
        "results/variants/core/core.full.aln"
    output:
        tree="results/phylogeny/raw.treefile",
        report="results/phylogeny/raw.iqtree",
    params:
        model=config["phylogeny"]["model"],
        bootstrap=config["phylogeny"]["bootstrap"],
        alrt=config["phylogeny"]["alrt"],
        seed=config["phylogeny"]["seed"],
    threads:
        config["resources"]["iqtree_threads"]
    resources:
        mem_mb=config["resources"].get("iqtree_mem_mb", 8000),
    log:
        "logs/phylogeny/raw_iqtree.log"
    conda:
        "../envs/phylogeny.yaml"
    shell:
        r"""
        mkdir -p results/phylogeny $(dirname {log})
        IQTREE=$(command -v iqtree2 || command -v iqtree)
        "$IQTREE" -s {input} -pre results/phylogeny/raw -m {params.model} \
          -B {params.bootstrap} --alrt {params.alrt} -T {threads} \
          -seed {params.seed} --redo > {log} 2>&1
        """
