rule integrated_report:
    input:
        REPORT_INPUTS
    output:
        html="results/report/index.html",
        summary="results/report/summary.json",
        config="results/report/config.resolved.json",
    params:
        samples=config["samples"],
        recombination=RECOMBINATION_ENABLED,
        accessory=ACCESSORY_ENABLED,
        functional=FUNCTIONAL_ENABLED,
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/report.py"
