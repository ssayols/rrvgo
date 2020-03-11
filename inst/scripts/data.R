#' prepareDemoData
#' 
#' Prepare demo data included in the package.
#' 
#' @return Whatever the call to write.table() would return
#' @details Taken asis from https://yulab-smu.github.io/clusterProfiler-book/chapter3.html#input-data,
#' which at the same time it is taken from the R package breastCancerMAINZ. It
#' contains 200 samples with breast cancer at different grades (I, II and III).
#' The dataset basically contains log2 ratios of the geometric means of
#' grade III vs. grade I samples.
#' @importFrom utils data write.table
prepareDemoData <- function() {
  data(geneList, package="DOSE")
  gene <- names(geneList)[abs(geneList) > 2]
  orgdb <- rrvgo:::loadOrgdb("org.Hs.eg.db")
  ego <- clusterProfiler::enrichGO(gene          = gene,
                                   universe      = names(geneList),
                                   OrgDb         = orgdb,
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
  write.table(subset(ego@result, qvalue < .05), file="inst/extdata/example.txt", sep="\t")
}
