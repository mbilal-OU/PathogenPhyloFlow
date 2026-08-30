# Contributing

Contributions that improve scientific validity, reproducibility, portability, testing, documentation, or visualization are welcome.

## Before opening a pull request

1. Keep analytical choices explicit in configuration.
2. Preserve raw outputs when a filtering step changes biological interpretation.
3. Add or update tests for new Python logic.
4. Document new output files and assumptions.
5. Avoid universal outbreak or transmission thresholds unless they are user supplied and clearly contextualized.
6. Run `pytest -q` and `snakemake --lint --snakefile workflow/Snakefile` when possible.

Bug reports should include the PathogenPhyloFlow commit, command, configuration, sample count, operating system, and relevant log excerpt.
