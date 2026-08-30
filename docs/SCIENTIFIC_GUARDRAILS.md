# Scientific guardrails

PathogenPhyloFlow is designed to automate computation without automating biological overclaiming.

## SNP distance is not transmission

Pairwise SNP distance can support relatedness assessment, but there is no universal SNP cutoff that proves direct transmission. Mutation rate, sampling interval, within-host diversity, recombination, organism, lineage, and surveillance design all matter. The workflow therefore reports candidate genomic proximity rather than declaring transmission chains.

## Recombination is measured, not assumed

Homologous recombination can create dense SNP blocks that distort clonal phylogeny. Its prevalence differs across taxa, lineages, and datasets. When recombination analysis is enabled, PathogenPhyloFlow keeps the raw result and the recombination-aware result so the effect can be inspected.

The current masking strategy is intentionally conservative: alignment intervals implicated by Gubbins are masked across the alignment before the recombination-aware IQ-TREE analysis. The masked fraction is reported.

## Reference choice can change the result

Reference-based variant calling can lose information in regions absent or highly divergent from the reference. Automatic reference assessment therefore uses a transparent medoid-style Mash criterion after basic assembly-quality screening. The selected reference and all candidate metrics are saved.

For datasets with strong population structure, users should inspect whether a single reference is appropriate.

## Temporal signal is screened before clock inference

Root-to-tip regression is useful as a diagnostic but is not sufficient evidence for a reliable molecular clock. PathogenPhyloFlow can prevent automatic TreeTime execution when the configured minimum sample count, sampling span, positive slope, and R-squared threshold are not met.

Publication-grade temporal inference may also require date randomization, sensitivity analyses, clock-model comparison, and epidemiological justification.

## Accessory content matters

Two isolates can be close in core SNP space but differ in plasmids, genomic islands, resistance determinants, virulence loci, and other accessory genes. The SNP/accessory discordance module flags pairs where low SNP distance co-occurs with unusually large accessory-gene distance.

A flagged pair is an interpretation target, not a biological conclusion.

## Functional hits require validation

AMR and virulence database matches depend on database content, identity thresholds, gene fragmentation, assembly quality, and nomenclature. Functional screening is therefore optional and raw hit tables are retained.

## Intended use

The workflow is intended for research, surveillance analysis, and reproducible exploratory genomics. It does not replace case investigation, clinical interpretation, infection-control decision making, or public-health policy.
