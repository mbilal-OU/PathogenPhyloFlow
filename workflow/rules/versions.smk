# Software-version provenance.
#
# Each rule below runs inside one of the workflow's Conda environments and
# records its tool's version string into results/versions/<tool>.txt. The
# collect_software_versions rule aggregates them into
# results/report/software_versions.json. Version capture never fails the
# workflow: an unresolvable version is recorded as "unknown".

VERSION_TOOLS = [
    "datasets",
    "mash",
    "snippy",
    "gubbins",
    "snp-dists",
    "iqtree",
    "prokka",
    "panaroo",
    "treetime",
]
if FUNCTIONAL_ENABLED:
    VERSION_TOOLS.append("abricate")


rule versions_input_env:
    output:
        "results/versions/datasets.txt",
    conda:
        "../envs/input.yaml"
    shell:
        "mkdir -p results/versions && (datasets --version > {output} 2>&1 || echo 'unknown' > {output})"


rule versions_reference_env:
    output:
        "results/versions/mash.txt",
    conda:
        "../envs/reference.yaml"
    shell:
        "mkdir -p results/versions && (mash --version > {output} 2>&1 || echo 'unknown' > {output})"


rule versions_snippy_env:
    output:
        "results/versions/snippy.txt",
    conda:
        "../envs/snippy.yaml"
    shell:
        "mkdir -p results/versions && (snippy --version > {output} 2>&1 || echo 'unknown' > {output})"


rule versions_recombination_env:
    output:
        "results/versions/gubbins.txt",
    conda:
        "../envs/recombination.yaml"
    shell:
        "mkdir -p results/versions && (run_gubbins.py --version > {output} 2>&1 || echo 'unknown' > {output})"


rule versions_phylogeny_env:
    output:
        dists="results/versions/snp-dists.txt",
        iqtree="results/versions/iqtree.txt",
    conda:
        "../envs/phylogeny.yaml"
    shell:
        r"""
        mkdir -p results/versions
        (snp-dists -v > {output.dists} 2>&1 || echo 'unknown' > {output.dists})
        IQTREE=$(command -v iqtree2 || command -v iqtree || echo iqtree2)
        ("$IQTREE" --version > {output.iqtree} 2>&1 || echo 'unknown' > {output.iqtree})
        """


rule versions_prokka_env:
    output:
        "results/versions/prokka.txt",
    conda:
        "../envs/prokka.yaml"
    shell:
        "mkdir -p results/versions && (prokka --version > {output} 2>&1 || echo 'unknown' > {output})"


rule versions_panaroo_env:
    output:
        "results/versions/panaroo.txt",
    conda:
        "../envs/panaroo.yaml"
    shell:
        "mkdir -p results/versions && (panaroo --version > {output} 2>&1 || echo 'unknown' > {output})"


rule versions_temporal_env:
    output:
        "results/versions/treetime.txt",
    conda:
        "../envs/temporal.yaml"
    shell:
        "mkdir -p results/versions && (treetime --version > {output} 2>&1 || echo 'unknown' > {output})"


if FUNCTIONAL_ENABLED:

    rule versions_functional_env:
        output:
            "results/versions/abricate.txt",
        conda:
            "../envs/functional.yaml"
        shell:
            "mkdir -p results/versions && (abricate --version > {output} 2>&1 || echo 'unknown' > {output})"


rule collect_software_versions:
    input:
        expand("results/versions/{tool}.txt", tool=VERSION_TOOLS),
    output:
        json="results/report/software_versions.json",
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/collect_versions.py"
