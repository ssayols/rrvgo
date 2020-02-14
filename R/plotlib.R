#' scatterPlot
#' Plot GO terms as scattered points.
#' 
#' @param simMatrix a (square) similarity matrix.
#' @param reducedTerms a data.frame with the reduced terms from reduceSimMatrix()
#' @param size what to use as point size. Can be either GO term's "size" or "score"
#' @param addLabel add labels with the most representative term of the group.
#' @param labelSize text size in the label.
#' @details  Distances between points represent the similarity between terms.
#' Axes are the first 2 components of applying a PCoA to the (di)similarity matrix.
#' Size of the point represents the provided scores or, in its absence, the number
#' of genes the GO term contains.
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
#' simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
#' reduced_go_analysis <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' scatterPlot(simMatrix, reduced_go_analysis)
#' @import ggplot2
#' @import ggrepel
#' @importFrom grid unit
#' @export
scatterPlot <- function(simMatrix, reducedTerms, size="score", addLabel=TRUE, labelSize=3) {

  if(!all(sapply(c("ggplot2", "ggrepel"), requireNamespace, quietly=TRUE))) {
    stop("Packages ggplot2, ggrepel and/or its dependencies not available. ",
         "Consider installing them before using this function.", call.=FALSE)
  }

  x <- cmdscale(as.matrix(as.dist(1-simMatrix)), eig=TRUE, k=2)

  df <- as.data.frame(x$points)
  df$term   <- reducedTerms$term[match(rownames(df), reducedTerms$go)]
  df$parent <- reducedTerms$parent[match(rownames(df), reducedTerms$go)]
  df$size   <- reducedTerms[match(rownames(df), reducedTerms$go), size]

  p <-
    ggplot2::ggplot(df, ggplot2::aes(x=V1, y=V2, label=term, color=parent)) +
      ggplot2::geom_point(ggplot2::aes(size=size), alpha=.5) +
      ggplot2::scale_color_discrete(guide=FALSE) +
      ggplot2::scale_size_continuous(guide=FALSE, range=c(0, 25)) +
      ggplot2::scale_x_continuous(name="") +
      ggplot2::scale_y_continuous(name="") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank())
  
  if(addLabel) {
    p + ggrepel::geom_label_repel(data=subset(df, parent == rownames(df)), box.padding=grid::unit(1, "lines"), size=labelSize)
  } else {
    p
  }
}

#' treemapPlot
#' Plot GO terms as a treemap.
#' 
#' @param reducedTerms a data.frame with the reduced terms from reduceSimMatrix()
#' @param size what to use as point size. Can be either GO term's "size" or "score"
#' @param ... other parameters sent to treemap::treemap()
#' @details  Distances between points represent the similarity between terms.
#' Axes are the first 2 components of applying a PCA to the similarity matrix.
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
#' simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
#' reduced_go_analysis <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' treemapPlot(reduced_go_analysis)
#' @importFrom treemap treemap
#' @export
treemapPlot <- function(reducedTerms, size="score", ...) {
  if(!all(sapply(c("treemap"), requireNamespace, quietly=TRUE))) {
    stop("Package treemap and/or its dependencies not available. ",
         "Consider installing it before using this function.", call.=FALSE)
  }

  reducedTerms$parent <- ifelse(reducedTerms$parent == "", reducedTerms$go, reducedTerms$parent)
  
  treemap::treemap(reducedTerms, index=c("parent", "term"), vSize=size, type="index", 
                   fontcolor.labels=c("#FFFFFFDD", "#00000080"), bg.labels=0, border.col="#00000080", ...)
}

#' wordlcoudPlot
#' Plot GO reduced terms as a wordcloud.
#' 
#' @param reducedTerms a data.frame with the reduced terms from reduceSimMatrix().
#' @param onlyParents use only parent terms to calculate frequencies.
#' @param ... other parameters sent to wordcloud::wordcloud()
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
#' simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
#' reduced_go_analysis <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' wordcloudPlot(reduced_go_analysis, min.freq=2, colors="black")
#' @importFrom tm Corpus TermDocumentMatrix
#' @importFrom wordcloud wordcloud
#' @export
wordcloudPlot <- function(reducedTerms, onlyParents=TRUE, ...) {
  if(!all(sapply(c("wordcloud", "tm", "slam"), requireNamespace, quietly=TRUE))) {
    stop("Package wordcloud and/or its dependencies (tm, slam) not available. ",
         "Consider installing it before using this function.", call.=FALSE)
  }
    
  if(onlyParents) {
    x <- tm::Corpus(tm::VectorSource(reducedTerms$term[reducedTerms$parent == ""]))
  } else {
    x <- tm::Corpus(tm::VectorSource(reducedTerms$term))
  }
  tdm <- tm::TermDocumentMatrix(x, control=list(removePunctuation=TRUE,
                                                stopwords=TRUE))
  m <- as.matrix(tdm)
  v <- sort(rowSums(m), decreasing=TRUE)
  d <- data.frame(word=names(v), freq=v)
  wordcloud::wordcloud(d$word, d$freq, ...)
}

#' heatmapPlot
#' Plot similarity matrix as a heatmap
#' 
#' @param simMatrix a (square) similarity matrix.
#' @param reducedTerms a data.frame with the reduced terms from reduceSimMatrix()
#' @param annotateParent whether to add annotation of the parent
#' @param annotationLabel display "go" ids or go "term" string
#' @param ... other parameters sent to pheatmap::pheatmap()
#' @details  Matrix with similarity scores between terms is represented as a heatmap.
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
#' simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
#' reduced_go_analysis <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' heatmapPlot(simMatrix, reduced_go_analysis, annotateParent=TRUE, annotationLabel="term", fontsize=6)
#' @importFrom pheatmap pheatmap
#' @export
heatmapPlot <- function(simMatrix, reducedTerms=NULL, annotateParent=FALSE, annotationLabel="go", ...) {
  
  if(!all(sapply(c("pheatmap"), requireNamespace, quietly=TRUE))) {
    stop("Package pheatmap and/or its dependencies not available. ",
         "Consider installing them before using this function.", call.=FALSE)
  }
  
  if(annotateParent && is.null(reducedTerms)) {
    warning("Need to provide a reducedTerms data.frame from reduceSimMatrix() to annotate the heatmap with parent terms.")
  }
  
  if(annotateParent && !is.null(reducedTerms)) {
    reducedTerms$ann <- reducedTerms$term[match(reducedTerms$parent, reducedTerms[, annotationLabel])]
    ann <- data.frame(parent=reducedTerms$ann[match(rownames(simMatrix), reducedTerms$go)],
                      row.names=rownames(simMatrix))
    pheatmap::pheatmap(simMatrix, annotation_row=ann, ...)
  } else {
    pheatmap::pheatmap(simMatrix, ...)
  }
  
}