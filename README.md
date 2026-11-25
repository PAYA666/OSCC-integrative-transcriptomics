# OSCC integrative transcriptomics

Source code and selected derived results accompanying **Integrative transcriptomic analysis identifies smoking and cisplatin resistance programs associated with clinical outcomes in oral squamous cell carcinoma**.

## Contents

- `historical/`: staged Part 1, Part 2 and Part 3 source recovered from the authors' project archives.
- `revision/`: analysis-specific revision scripts, including strict nested cross-validation, smoking sensitivities, independent cohort evaluation, pathway checks and figure generation.
- `parameters/TCGA_reference.tsv`: frozen Final8 genes, coefficients, centres and standard deviations.
- `results/`: gene-level, technical-reference and aggregate outputs associated with the revised manuscript. These are selected tables, not the full journal supplement.
- `docs/`: run order, input requirements, file provenance and exact-duplicate table aliases.

## Start here

Read [RUNNING.md](docs/RUNNING.md) before executing scripts and [DATA_SOURCES.md](docs/DATA_SOURCES.md) to acquire the required public source data. Files are organized by analysis, not as one monolithic command. Use separate R sessions and a writable analysis workspace: some historical scripts reset the R environment or replace their own derived outputs.

The revised prognostic analysis is in `revision/02_strict_cv/`; historical apparent/fixed-signature estimates are not the primary internally validated performance estimates. The strict pipeline repeats feature selection within outer training folds. External frozen Final8 scores use cohort-specific expression processing as described in the manuscript.

## Scope and reproducibility

This release provides staged source and selected frozen results. Original matrices, patient-level prediction objects, large intermediate RDS files, and separately distributed annotation/gene-set collections are not bundled. Historical scripts may require intermediate objects generated in earlier project stages; see the run-order and input notes. A clean-environment end-to-end execution of this assembled release has not been completed. Python syntax and R lexical structure were checked during packaging; these checks are not runtime validation.

The Part 3 convergence reconstruction is retained as a reconstructed source candidate, not relabeled as a recovered original. Source provenance is indexed in `docs/file_provenance.tsv`. Code comments and docstrings were removed at the authors' request; operating instructions are retained in this documentation. Workstation data-root literals were replaced by environment-variable references without intentionally changing analytical formulas.

## Interpretation

Final8 is a research prognostic score, not a validated single-specimen clinical assay. SEQC results concern reference RNA mixtures, not paired OSCC patients. GSE84846 results evaluate deposited binary codes, not censoring-based outcomes. E-MTAB-15148 summaries describe one SCC-9 selection series and do not treat cells as independent biological replicates.

## Licensing and citation

No open-source license has been selected for this snapshot. Public availability alone does not grant an additional software reuse license. Third-party data and resources remain subject to their original terms. A code license can be added after author confirmation.

The associated manuscript is under revision. No accepted-publication DOI or release DOI is claimed here. When referring to this repository, record the repository URL and the exact commit used.
