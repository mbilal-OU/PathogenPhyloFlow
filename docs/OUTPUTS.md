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
```

These paths are intended to remain stable across minor releases whenever possible.
