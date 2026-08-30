if RECOMBINATION_ENABLED:

    rule gubbins:
        input:
            "results/variants/core/core.full.aln"
        output:
            gff="results/recombination/gubbins.recombination_predictions.gff",
            filtered="results/recombination/gubbins.filtered_polymorphic_sites.fasta",
            tree="results/recombination/gubbins.final_tree.tre",
        threads:
            config["resources"]["gubbins_threads"]
        log:
            "logs/recombination/gubbins.log"
        conda:
            "../envs/recombination.yaml"
        shell:
            r"""
            mkdir -p results/recombination $(dirname {log})
            run_gubbins.py --prefix results/recombination/gubbins \
              --threads {threads} {input} > {log} 2>&1
            """


    rule mask_recombination:
        input:
            alignment="results/variants/core/core.full.aln",
            gff="results/recombination/gubbins.recombination_predictions.gff",
        output:
            alignment="results/recombination/core.masked.full.aln",
            summary="results/recombination/mask_summary.json",
        conda:
            "../envs/recombination.yaml"
        script:
            "../scripts/mask_recombination.py"


    rule filtered_snp_distances:
        input:
            "results/recombination/core.masked.full.aln"
        output:
            "results/recombination/filtered.snp_dist.tsv"
        conda:
            "../envs/phylogeny.yaml"
        shell:
            "snp-dists {input} > {output}"


    rule filtered_iqtree:
        input:
            "results/recombination/core.masked.full.aln"
        output:
            tree="results/phylogeny/recombination_filtered.treefile",
            report="results/phylogeny/recombination_filtered.iqtree",
        params:
            model=config["phylogeny"]["model"],
            bootstrap=config["phylogeny"]["bootstrap"],
            alrt=config["phylogeny"]["alrt"],
            seed=config["phylogeny"]["seed"],
        threads:
            config["resources"]["iqtree_threads"]
        log:
            "logs/phylogeny/recombination_filtered_iqtree.log"
        conda:
            "../envs/phylogeny.yaml"
        shell:
            r"""
            mkdir -p results/phylogeny $(dirname {log})
            IQTREE=$(command -v iqtree2 || command -v iqtree)
            "$IQTREE" -s {input} -pre results/phylogeny/recombination_filtered \
              -m {params.model} -B {params.bootstrap} --alrt {params.alrt} \
              -T {threads} --seed {params.seed} --redo > {log} 2>&1
            """


    rule compare_trees:
        input:
            raw="results/phylogeny/raw.treefile",
            filtered="results/phylogeny/recombination_filtered.treefile",
        output:
            json="results/phylogeny/tree_comparison.json",
            tsv="results/phylogeny/tree_comparison.tsv",
        conda:
            "../envs/report.yaml"
        script:
            "../scripts/tree_compare.py"


    rule finalize_phylogeny:
        input:
            tree="results/phylogeny/recombination_filtered.treefile",
            alignment="results/recombination/core.masked.full.aln",
            distances="results/recombination/filtered.snp_dist.tsv",
        output:
            tree="results/phylogeny/final.treefile",
            alignment="results/phylogeny/final.aln",
            distances="results/phylogeny/final.snp_dist.tsv",
        shell:
            "cp {input.tree} {output.tree} && cp {input.alignment} {output.alignment} && cp {input.distances} {output.distances}"

else:

    rule finalize_phylogeny:
        input:
            tree="results/phylogeny/raw.treefile",
            alignment="results/variants/core/core.full.aln",
            distances="results/variants/core/raw.snp_dist.tsv",
        output:
            tree="results/phylogeny/final.treefile",
            alignment="results/phylogeny/final.aln",
            distances="results/phylogeny/final.snp_dist.tsv",
        shell:
            "cp {input.tree} {output.tree} && cp {input.alignment} {output.alignment} && cp {input.distances} {output.distances}"
