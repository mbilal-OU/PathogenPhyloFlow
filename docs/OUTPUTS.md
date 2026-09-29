# Output reference

PathogenPhyloFlow writes each analytical layer to a separate directory so provenance remains clear.

| Directory | Main contents | Interpretation |
|---|---|---|
| `results/input/` | Resolved assemblies and manifest records | Exact sequence input used by the run |
| `results/reference/` | Selected reference, candidate metrics, selection JSON | Why a reference was selected |
| `results/variants/` | Per-isolate Snippy outputs, core alignment, SNP matrices | Reference-based core variation |
| `results/recombination/` | Gubbins outputs, masked alignment, masking summary | Recombination evidence and masking impact |
| `results/phylogeny/` | Raw, recombination-aware, and final trees | Evolutionary relationships before and after filtering |
| `results/accessory/` | Prokka annotations and Panaroo matrices | Accessory-gene distribution |
| `results/functional/` | Optional ABRicate results | Candidate AMR and virulence determinants |
| `results/integration/` | SNP/accessory discordance tables | Isolate pairs with contrasting core and accessory similarity |
| `results/temporal/` | Root-to-tip metrics and optional TreeTime output | Temporal screening and clock analysis status |
| `results/report/` | HTML report and machine-readable summary | Integrated run overview |
| `results/versions/` | Per-tool version strings | Software provenance for the run |

## Stable high-level files

The workflow exposes a small number of stable files for downstream use:

```text
results/phylogeny/final.treefile
results/phylogeny/final.aln
results/phylogeny/final.snp_dist.tsv
results/integration/snp_accessory_discordance.tsv
results/temporal/temporal_screen.json
results/report/index.html
results/report/summary.json
results/report/software_versions.json
```

These paths are intended to remain stable across minor releases whenever possible.

### A note on the `Reference` tip

`snippy-core` includes the reference genome itself as a taxon, so
`final.treefile` (and the SNP distance matrices) contain N samples **plus one
`Reference` tip** paired with whichever genome was auto-selected as the
reference. Treat it as an explicit outgroup-like taxon, not as an additional
sample: sample counts in the report refer to the N input samples, and the
temporal screen silently skips the undated `Reference` tip.

### Software provenance

`results/report/software_versions.json` records the version string of every
external tool (one per Conda environment) plus the Python, Snakemake, and git
commit the run used. The raw per-tool files live in `results/versions/`.
When reporting results, cite these versions alongside the workflow version.
