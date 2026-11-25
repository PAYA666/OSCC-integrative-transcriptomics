# OSCC integrative transcriptomics

Analysis code and selected derived results for **Integrative transcriptomic analysis identifies smoking and cisplatin resistance programs associated with clinical outcomes in oral squamous cell carcinoma**.

## Repository contents

| Directory | Contents |
|---|---|
| `historical/` | Archived Part 1, Part 2 and Part 3 analysis scripts |
| `revision/` | Strict nested cross-validation, sensitivity analyses, external evaluation and figure scripts |
| `parameters/` | Frozen Final8 coefficients and reference scaling parameters; published comparator coefficients |
| `results/` | Selected gene-level, technical-reference and aggregate results |
| `docs/` | Data sources, execution guidance, analysis mapping and file provenance |

Start with [execution guidance](docs/RUNNING.md), [data sources](docs/DATA_SOURCES.md) and the [analysis map](docs/ANALYSIS_MAP.md). The scripts are organized by analysis stage and should be run in a separate writable workspace with their required inputs. Historical stages are listed in [run order](docs/historical_run_order.tsv).

## Analyses and interpretation

The revised internal prognostic analysis is in `revision/02_strict_cv/`. It repeats feature selection within outer training folds; historical apparent and fixed-signature estimates are not the primary internally validated estimates. External Final8 evaluation uses frozen coefficients with the cohort-specific expression processing described in the manuscript.

Final8 is a research prognostic score, not a validated single-specimen clinical assay. SEQC evaluates reference RNA mixtures rather than paired OSCC patients. GSE84846 evaluates deposited binary outcome codes rather than survival time and censoring. E-MTAB-15148 provides descriptive expression trajectories from one SCC-9 selection series, without treating cells as independent biological replicates.

## Inputs and reproducibility

Original expression matrices, patient-level prediction objects, large intermediate RDS files and separately distributed annotation or gene-set collections are not included. Obtain the source data under their applicable terms and supply the intermediate objects required by each stage. The `results/` directory contains selected tables, not a copy of the complete journal supplement; duplicate-table aliases are documented in [result_aliases.tsv](docs/result_aliases.tsv).

Archived scripts include a reconstructed Part 3 convergence stage; consult the provenance records when using it. Workstation-specific data roots are configured through environment variables. Operating instructions and validation records are in `docs/`.

The assembled release has not undergone a complete clean-environment rerun. Packaging checks comprise Python syntax checks and R lexical delimiter checks. The four-reference SEQC weighted-profile script was additionally executed with the included input tables and reproduced the frozen S30 result byte-for-byte. See [RUNNING.md](docs/RUNNING.md) for the scope of these checks.

## License and citation

No open-source software license has yet been selected for this repository; third-party resources retain their original terms. The associated manuscript is under revision. Cite the repository URL and the exact commit used when referring to this version.
