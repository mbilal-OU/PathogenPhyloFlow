rule temporal_analysis:
    input:
        tree="results/phylogeny/final.treefile",
        alignment="results/phylogeny/final.aln",
        samples=config["samples"],
    output:
        screen="results/temporal/temporal_screen.json",
        points="results/temporal/root_to_tip.tsv",
        status="results/temporal/treetime_status.txt",
    params:
        settings=config["temporal"]
    log:
        "logs/temporal/temporal.log"
    conda:
        "../envs/temporal.yaml"
    script:
        "../scripts/temporal_analysis.py"
