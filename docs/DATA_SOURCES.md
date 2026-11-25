# Data sources

| Context | Source |
|---|---|
| Tumor, smoking, model development | [TCGA-HNSC through GDC](https://portal.gdc.cancer.gov/projects/TCGA-HNSC) |
| Resistant cell models | GEO GSE197561, GSE111585 |
| External survival evaluation | GEO GSE41613, GSE42743, GSE85446 |
| Tumor single-cell context | GEO GSE103322 |
| Exploratory binary-status evaluation | GEO GSE84846 |
| Reference RNA assay comparison | GEO GSE47774 and GSE56457; GPL10558 |
| SCC-9 selection-series expression | [E-MTAB-15148](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15148) |

GEO accessions can be resolved at `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=ACCESSION`.

Obtain original matrices, annotations and clinical data from the source repositories under their applicable terms. Acquire the MSigDB Hallmark collection separately; the analysis used h.all.v2026.1.Hs.symbols.gmt. No complete gene-set library is redistributed here. Summarized enrichment results are derived analysis outputs.

The downloadable results in this repository omit patient-level prediction files, donor-prefix score tables and combined patient-level RDS objects. Cell-line GSM labels in leave-one-replicate-out outputs identify public experimental replicates, not patients. Numeric result files are unchanged from the selected source tables; duplicate byte-identical files are mapped in result_aliases.tsv.

## Input filenames expected by the SEQC script

- `GPL10558_SOFT.txt` (platform annotation)
- `GSE56457_SEQC_qPCR_GEOSub.txt.gz` (PrimePCR data)
- `GSE47774_SEQC_ILM_COH.txt.gz` (RNA-seq counts)
- `GSE56457_SEQC_ILM_GEOSub.txt.gz` (HT-12 array data)

Place these four files in the directory passed as `--input-dir`. They are not part of this release. This source-data list is not an instruction to download files automatically.
