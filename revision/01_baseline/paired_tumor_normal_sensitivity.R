suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

root <- Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA"))
meta <- fread("results/sample_manifest.tsv", na.strings = "NA")
normal <- meta[sample_type == "Solid Tissue Normal", unique(patient_barcode)]
paired <- meta[patient_barcode %in% normal &
                 sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
counts_meta <- paired[, .(n_tumor = sum(sample_type == "Primary Tumor"),
                          n_normal = sum(sample_type == "Solid Tissue Normal")),
                      by = patient_barcode]
stopifnot(nrow(counts_meta) == 32L, all(counts_meta$n_tumor == 1L),
          all(counts_meta$n_normal == 1L), nrow(paired) == 64L)

raw <- readRDS(file.path(root, "TCGA_OSCC_Nature8_raw_count_matrix.rds"))
raw <- raw[raw$gene_type == "protein_coding", ]
cnt <- as.matrix(raw[, paired$sample_barcode])
rownames(cnt) <- raw$gene_id
storage.mode(cnt) <- "integer"
coldata <- data.frame(
  patient = factor(paired$patient_barcode),
  type = factor(ifelse(paired$sample_type == "Primary Tumor", "Tumor", "Normal"),
                levels = c("Normal", "Tumor")),
  row.names = paired$sample_barcode
)
stopifnot(identical(colnames(cnt), rownames(coldata)))
keep <- rowSums(cnt >= 10L) >= 10L
dds <- DESeqDataSetFromMatrix(cnt[keep, , drop = FALSE], coldata,
                              design = ~ patient + type)
dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)
paired_result <- as.data.table(results(dds, contrast = c("type", "Tumor", "Normal"),
                                      alpha = 0.05), keep.rownames = "gene_id")

original <- fread(file.path(root, "TCGA_OSCC_Nature8_DESeq2_ALL.tsv"),
                  na.strings = "NA")
stopifnot(!anyDuplicated(original$gene_id),
          sum(is.finite(original$padj) & original$padj < .05 &
              abs(original$log2FoldChange) > 1, na.rm=TRUE) == 4461L)
comparison <- merge(
  paired_result[, .(gene_id, paired_log2FC = log2FoldChange,
                    paired_padj = padj)],
  original[, .(gene_id, unpaired_log2FC = log2FoldChange,
               unpaired_padj = padj)],
  by = "gene_id", all = FALSE
)
comparison[, `:=`(
  paired_strict = is.finite(paired_padj) & paired_padj < 0.05 &
    abs(paired_log2FC) > 1,
  unpaired_strict = is.finite(unpaired_padj) & unpaired_padj < 0.05 &
    abs(unpaired_log2FC) > 1,
  same_direction = sign(paired_log2FC) == sign(unpaired_log2FC)
)]

summary <- data.table(
  paired_patients = nrow(counts_meta), paired_samples = nrow(paired),
  paired_strict_DEGs = sum(comparison$paired_strict, na.rm = TRUE),
  original_unpaired_strict_DEGs = sum(comparison$unpaired_strict, na.rm = TRUE),
  overlapping_strict_DEGs = sum(comparison$paired_strict &
                                  comparison$unpaired_strict, na.rm = TRUE),
  same_direction_among_overlapping = sum(comparison$paired_strict &
    comparison$unpaired_strict & comparison$same_direction, na.rm = TRUE),
  spearman_log2FC = cor(comparison$paired_log2FC, comparison$unpaired_log2FC,
                        method = "spearman", use = "complete.obs"),
  pearson_log2FC = cor(comparison$paired_log2FC, comparison$unpaired_log2FC,
                       method = "pearson", use = "complete.obs")
)

fwrite(paired_result, "results/paired_tumor_normal_DESeq2.tsv", sep = "\t")
fwrite(comparison, "results/paired_tumor_normal_comparison.tsv", sep = "\t")
fwrite(summary, "results/paired_tumor_normal_summary.tsv", sep = "\t")
capture.output(sessionInfo(), file = "results/paired_tumor_normal_session.txt")
