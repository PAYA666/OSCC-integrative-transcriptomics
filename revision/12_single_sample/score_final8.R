score_final8 <- function(x, reference, expression_scale) {
  if (!identical(expression_scale, "TCGA_Nature8_VST"))
    stop("Input must already be on the archived TCGA Nature8 VST scale.")
  stopifnot(nrow(reference)==8L, !anyDuplicated(reference$gene),
            length(x)==8L, is.numeric(x), all(is.finite(x)),
            !anyDuplicated(names(x)), setequal(names(x), reference$gene),
            all(is.finite(as.matrix(reference[c("center","sd","beta")]))),
            all(reference$sd>0))
  sum((x[reference$gene]-reference$center)/reference$sd * reference$beta)
}
