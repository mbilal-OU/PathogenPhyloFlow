# Environment setup

PathogenPhyloFlow uses `environment.yaml` for the workflow-level environment.

## Preferred: Mamba

If Mamba is already available in your base Conda installation:

```bash
mamba env create -f environment.yaml
conda activate pathogenphyloflow
```

## Conda-only fallback

Mamba is not required to bootstrap the project. If `mamba` is not installed yet, create the same environment directly with Conda:

```bash
conda env create -f environment.yaml
conda activate pathogenphyloflow
```

After activation, confirm that Snakemake is coming from the project environment before launching a workflow:

```bash
snakemake --version
```

The workflow itself should still be run with `--use-conda` so rule-specific environments are resolved as defined in `workflow/envs/`.
