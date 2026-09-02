# Sample sheet requirements

PathogenPhyloFlow reads the sample sheet as a tab-delimited TSV file. The workflow validates the sheet before constructing the analysis DAG, so malformed inputs fail early.

## Required columns

The following columns must be present:

```text
sample	assembly	accession	date	location	host
```

Additional metadata columns are allowed and are retained for reporting.

## Input rules

- The sheet must contain at least two sample rows.
- `sample` identifiers must be unique.
- `sample` identifiers may contain only letters, numbers, underscores, periods, and hyphens.
- Each row must provide exactly one genome source: either a local path in `assembly` or an NCBI assembly accession in `accession`.
- Do not populate both `assembly` and `accession` for the same sample, and do not leave both blank.
- Keep the file tab-delimited; comma-separated input is not parsed as a valid sample sheet.

A minimal mixed-input example is available at [`config/samples.example.tsv`](../config/samples.example.tsv).

## Validate before a full run

A dry run exercises sample-sheet loading and validation without launching the full analysis:

```bash
snakemake --snakefile workflow/Snakefile --use-conda --cores 1 --dry-run
```

If validation fails, correct the reported sample identifier, missing column, duplicate row, or genome-source field before starting the workflow.
