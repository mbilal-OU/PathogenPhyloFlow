# E. coli O157:H7 real-data validation set

This example is a small, curated set of 10 public *Escherichia coli* O157:H7 assemblies selected to exercise the major PathogenPhyloFlow branches without requiring a large surveillance dataset.

It intentionally mixes well-known human/food outbreak isolates with cattle-associated isolates. That gives the workflow enough biological contrast to test reference assessment, core-SNP phylogeny, recombination analysis, accessory-genome differences, and metadata-aware interpretation.

## Included isolates

| Sample | Assembly | Year used | Source/context |
|---|---|---:|---|
| EDL933 | GCF_000006665.1 | 1982 | Ground beef; hamburger-associated outbreak reference strain |
| Sakai | GCF_000008865.2 | 1996 | Human; Sakai outbreak reference strain |
| EC4115 | GCF_000021125.1 | 2006 | Human; spinach-associated outbreak |
| TW14359 | GCF_000022225.1 | 2006 | Human; spinach-associated outbreak |
| TW14588 | GCF_000155125.1 | 2006 | Lettuce-associated outbreak |
| Xuzhou21 | GCF_000262125.1 | 1999 | HUS patient; Xuzhou outbreak |
| SS17 | GCF_000730345.1 | not forced | Super-shedder cattle isolate |
| SS52 | GCF_000803705.1 | not forced | Super-shedder cattle isolate |
| JEONG-1266 | GCF_001558995.2 | not forced | Super-shedder cattle isolate |
| FRIK2069 | GCF_001651925.1 | 2011 | Cattle-associated isolate |

Dates are intentionally left blank when the literature gives an interval, an uncertain isolation year, or conflicting metadata. PathogenPhyloFlow should not invent precision merely to make a temporal analysis run.

## Why this is a validation dataset, not an outbreak reconstruction

These genomes do **not** represent one epidemiologically linked outbreak. They span different outbreaks, countries, sources, and years. The example therefore validates pipeline behavior and biological guardrails; it must not be interpreted as evidence of transmission among these isolates.

The expected high-level behavior is:

1. all 10 NCBI assembly accessions resolve reproducibly;
2. automatic reference assessment records a transparent recommendation;
3. Snippy generates a core alignment and pairwise SNP distances;
4. the raw SNP tree is retained;
5. Gubbins identifies and masks recombinant segments when supported by the data;
6. a second phylogeny is built from the recombination-aware alignment;
7. Panaroo measures accessory-genome variation;
8. SNP/accessory discordance is evaluated without treating a threshold as proof of transmission;
9. the temporal module reports a diagnostic only and does not claim clock-like evolution from root-to-tip regression alone;
10. the final HTML/JSON report keeps those evidence layers separate.

## Run it

Create the main environment:

```bash
mamba env create -f environment.yaml
conda activate pathogenphyloflow
```

Inspect the DAG:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile examples/ecoli_o157/config.yaml \
  --use-conda \
  --cores 8 \
  --dry-run
```

Run the analysis:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile examples/ecoli_o157/config.yaml \
  --use-conda \
  --cores 8
```

The main report is written to `results/report/index.html`.

Run this example in a clean working directory, or remove/relocate an existing `results/` directory first. PathogenPhyloFlow currently uses a common `results/` root for each run.

## Functional screening

The example keeps `functional.enabled: false` by default because ABRicate database installation/update state is external to the workflow and should not be silently assumed. After installing and verifying the desired databases, set:

```yaml
functional:
  enabled: true
  databases:
    - resfinder
    - vfdb
```

The functional branch can then be evaluated on the same genomes.

## Provenance used for curation

The accession/source choices were cross-checked against NCBI records and peer-reviewed genome/comparative studies, including:

- Perna NT et al. 2001. Genome sequence of enterohaemorrhagic *Escherichia coli* O157:H7. *Nature*.
- Hayashi T et al. 2001. Complete genome sequence of enterohemorrhagic *E. coli* O157:H7 Sakai. *DNA Research*.
- Kulasekara BR et al. 2009. Analysis of the genome of the 2006 spinach-associated outbreak isolate TW14359. *Infection and Immunity*.
- Eppinger M et al. 2011. Genome-scale analysis of the 2006 spinach-associated outbreak lineage including EC4115.
- Xiong Y et al. 2012. A novel *E. coli* O157:H7 clone causing a major HUS outbreak in China. *PLoS ONE*.
- Cote R et al. 2015 / related complete-genome reports for super-shedder isolates SS17 and SS52.
- Teng L et al. 2016. Complete genome sequence of JEONG-1266 from a super-shedder steer. *Genome Announcements*.
- Wang J et al. 2020. Core and accessory genome comparison of Australian and international O157 STEC. *Frontiers in Microbiology* 11:566415.
- Peroutka-Bigus N et al. 2024. Phenotypic and genomic comparison of outbreak and cattle-associated O157:H7 isolates. *Microbiology Spectrum* 12:e04140-23.

This file records the rationale for the tutorial set. NCBI accessions remain the authoritative sequence identifiers used by the workflow.
