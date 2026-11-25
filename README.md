# OSCC integrative transcriptomics

This repository accompanies our study, **Integrative transcriptomic analysis identifies smoking and cisplatin resistance programs associated with clinical outcomes in oral squamous cell carcinoma**. It provides analysis scripts, model parameters and selected derived results supporting the manuscript.

## Study objective and research contributions

We investigate how transcriptional programs associated with oral squamous cell carcinoma, smoking history and acquired cisplatin resistance converge or differ, and how selected expression features relate to clinical outcomes. Our approach combines comparisons across research contexts with co-expression analysis, single-cell expression profiling and prognostic evaluation.

The study brings together three complementary components:

1. **Tumor biology and prognosis:** characterization of tumor-associated expression and co-expression structure, followed by development and evaluation of the Final8 prognostic score.
2. **Smoking-associated expression:** comparison of smoking-history groups and assessment of available clinical, anatomical and cellular-composition covariates.
3. **Acquired cisplatin resistance:** comparison of independent resistant-cell models, prioritization of candidate genes and examination of their expression in additional cellular contexts.

The authors take responsibility for the study design, methodological choices, verification of reported results, and biological and clinical interpretation. Individual author contributions are described in the manuscript. This repository documents the computational methods and selected outputs so that readers can inspect the analyses underlying those scientific conclusions.

## Repository contents

| Directory | Contents |
|---|---|
| `historical/` | Archived Part 1, Part 2 and Part 3 analysis scripts |
| `revision/` | Strict nested cross-validation, sensitivity analyses, external evaluation and figure scripts |
| `parameters/` | Frozen Final8 coefficients and reference scaling parameters; published comparator coefficients |
| `results/` | Selected gene-level, technical-reference and aggregate results |
| `docs/` | Data sources, execution guidance, analysis mapping and file provenance |

Start with the [analysis map](docs/ANALYSIS_MAP.md), [data sources](docs/DATA_SOURCES.md) and [execution guidance](docs/RUNNING.md). Scripts follow the study's analytical stages; the original three-part workflow is indexed in [run order](docs/historical_run_order.tsv).

## Analyses and interpretation

The revised internal prognostic analysis is in `revision/02_strict_cv/`. We distinguish evaluation of the full model-building pipeline from evaluation of the frozen Final8 score: the internal pipeline repeats feature selection within outer training folds, while external evaluations use frozen coefficients with the cohort-specific expression processing described in the manuscript. Historical apparent and fixed-signature estimates are retained in the source archive rather than used as the primary internally validated estimates.

Each additional dataset serves a defined purpose. SEQC assesses technical correspondence across reference RNA assays; GSE84846 provides an exploratory evaluation of deposited binary outcome codes; and E-MTAB-15148 supplies descriptive expression trajectories from one SCC-9 selection series. These analyses complement, rather than replace, the external survival evaluations. Final8 is presented as a research prognostic score, and cells from a single selection series are not treated as independent biological replicates.

## Inputs and reproducibility

Run the scripts in a separate writable workspace with the inputs specified for each stage. Original expression matrices, patient-level prediction objects, large intermediate RDS files and separately distributed annotation or gene-set collections are not bundled. Obtain the source data under their applicable terms and supply the required intermediate objects. The `results/` directory contains selected tables rather than the complete journal supplement; duplicate-table aliases are documented in [result_aliases.tsv](docs/result_aliases.tsv).

Archived scripts include a reconstructed Part 3 convergence stage, identified in the provenance records. Data roots are configured through environment variables. Operating instructions and validation records are in `docs/`.

Release checks comprise Python syntax checks, R lexical delimiter checks and a local execution of the four-reference SEQC weighted-profile script, whose output reproduced the frozen S30 result byte-for-byte. The assembled release has not undergone a complete clean-environment rerun. See [RUNNING.md](docs/RUNNING.md) for input dependencies and the scope of these checks.

## License and citation

No open-source software license has yet been selected for this repository; third-party resources retain their original terms. The associated manuscript is under revision. Cite the repository URL and the exact commit used when referring to this version.
