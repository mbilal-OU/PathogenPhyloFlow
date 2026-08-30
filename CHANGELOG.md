# Changelog

All notable changes to PathogenPhyloFlow are documented here.

## 0.1.0 - 2026-08-30

Initial public release.

### Workflow

- Modular Snakemake architecture with separate rules for input handling, reference assessment, core-SNP analysis, recombination, accessory-genome analysis, temporal diagnostics, and reporting.
- Local assembly and NCBI assembly accession input support.
- Transparent reference assessment using assembly statistics and Mash distance.
- Core-SNP analysis with Snippy and phylogeny with IQ-TREE.
- Recombination analysis with Gubbins while retaining both raw and recombination-aware outputs.
- Pairwise SNP distance output.
- Prokka and Panaroo accessory-genome branch.
- Optional ABRicate-based ResFinder and VFDB screening.
- SNP versus accessory-genome discordance reporting.
- Root-to-tip temporal screening before optional TreeTime execution.
- Integrated HTML and JSON reporting.

### Reproducibility and engineering

- Installable `pathogenphyloflow` Python package for shared utilities.
- Python 3.11 pinned for the workflow runtime.
- Per-rule Conda environments.
- Configuration schema validation.
- Unit tests and Snakemake DAG validation.
- End-to-end four-genome workflow validation in GitHub Actions.
- Public NCBI accession-resolution smoke test.
- Deterministic handling of shared package imports across Snakemake rule environments.

### Documentation and examples

- Curated 10-genome *Escherichia coli* O157:H7 tutorial dataset.
- Four-genome CI validation subset with committed result preview.
- README figures generated from actual successful workflow outputs.
- Raw versus recombination-filtered tree comparison.
- SNP versus accessory-genome comparison figure.
- Temporal-screen diagnostic figure.
- Scientific guardrails and output reference documentation.
- `CITATION.cff`, MIT license, and contribution guidance.

### Interpretation scope

PathogenPhyloFlow reports genomic evidence layers without treating SNP proximity as proof of direct transmission. The bundled O157:H7 example is intended for workflow validation and teaching, not as a reconstruction of one epidemiologically linked outbreak.
