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
#' reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' scatterPlot(simMatrix, reducedTerms)
#' @import ggplot2
#' @import ggrepel
#' @importFrom grid unit
#' @importFrom stats as.dist cmdscale
#' @export
scatterPlot <- function(simMatrix, reducedTerms, size="score", addLabel=TRUE, labelSize=3) {

  if(!all(sapply(c("ggplot2", "ggrepel"), requireNamespace, quietly=TRUE))) {
    stop("Packages ggplot2, ggrepel and/or its dependencies not available. ",
         "Consider installing them before using this function.", call.=FALSE)
  }

  x <- cmdscale(as.matrix(as.dist(1-simMatrix)), eig=TRUE, k=2)

  df <- cbind(as.data.frame(x$points),
              reducedTerms[match(rownames(x$points), reducedTerms$go), c("term", "parent", "parentTerm", "size")])
  
  p <-
    ggplot2::ggplot(df, ggplot2::aes(x=V1, y=V2, color=parentTerm)) +
      ggplot2::geom_point(ggplot2::aes(size=size), alpha=.5) +
      ggplot2::scale_color_discrete(guide=FALSE) +
      ggplot2::scale_size_continuous(guide=FALSE, range=c(0, 25)) +
      ggplot2::scale_x_continuous(name="") +
      ggplot2::scale_y_continuous(name="") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank())
  
  if(addLabel) {
    p + ggrepel::geom_label_repel(aes(label=parentTerm),
                                  data=subset(df, parent == rownames(df)),
                                  box.padding=grid::unit(1, "lines"), size=labelSize)
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
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
#' simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
#' reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' treemapPlot(reducedTerms)
#' @importFrom treemap treemap
#' @export
treemapPlot <- function(reducedTerms, size="score", ...) {
  if(!all(sapply(c("treemap"), requireNamespace, quietly=TRUE))) {
    stop("Package treemap and/or its dependencies not available. ",
         "Consider installing it before using this function.", call.=FALSE)
  }

  treemap::treemap(reducedTerms, index=c("parentTerm", "term"), vSize=size, type="index", title="", 
                   palette=gg_color_hue(length(unique(reducedTerms$parent))),
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
#' reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' wordcloudPlot(reducedTerms, min.freq=1, colors="black")
#' @importFrom tm Corpus TermDocumentMatrix
#' @importFrom wordcloud wordcloud
#' @export
wordcloudPlot <- function(reducedTerms, onlyParents=TRUE, ...) {
  if(!all(sapply(c("wordcloud", "tm", "slam"), requireNamespace, quietly=TRUE))) {
    stop("Package wordcloud and/or its dependencies (tm, slam) not available. ",
         "Consider installing it before using this function.", call.=FALSE)
  }
    
  if(onlyParents) {
    x <- tm::Corpus(tm::VectorSource(reducedTerms$term[reducedTerms$parent == reducedTerms$go]))
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
#' @param annotationLabel display "parent" ids or "parentTerm" string
#' @param ... other parameters sent to pheatmap::pheatmap()
#' @details  Matrix with similarity scores between terms is represented as a heatmap.
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
#' simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
#' reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' heatmapPlot(simMatrix, reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=6)
#' @importFrom pheatmap pheatmap
#' @export
heatmapPlot <- function(simMatrix, reducedTerms=NULL, annotateParent=TRUE, annotationLabel="parentTerm", ...) {
  
  if(annotateParent && !annotationLabel %in% c("parent", "parentTerm")) {
    stop("annotationLabel should be one of c('parent', 'parentTerm')")
  }
  
  if(!all(sapply(c("pheatmap"), requireNamespace, quietly=TRUE))) {
    stop("Package pheatmap and/or its dependencies not available. ",
         "Consider installing them before using this function.", call.=FALSE)
  }
  
  if(annotateParent && is.null(reducedTerms)) {
    warning("Need to provide a reducedTerms data.frame from reduceSimMatrix() to annotate the heatmap with parent terms.")
  }
  
  if(annotateParent && !is.null(reducedTerms)) {
    ann <- data.frame(parent=factor(reducedTerms[match(rownames(simMatrix), reducedTerms$go), annotationLabel]),
                      row.names=rownames(simMatrix))
    annColors <- list(parent=gg_color_hue(length(unique(ann$parent))))
    names(annColors$parent) <- levels(ann$parent)
    pheatmap::pheatmap(simMatrix, annotation_row=ann, annotation_colors=annColors, ...)
  } else {
    pheatmap::pheatmap(simMatrix, ...)
  }
  
}

#' gg_color_hue
#' Emulate ggplot2 color palette.
#' 
#' @param n number of colors
#' @details  It is just equally spaced hues around the color wheel, starting from 15:
#' @importFrom grDevices hcl
#' @examples
#' \dontrun{
#' plot(1:10, pch=16, cex=2, col=gg_color_hue(10))
#' }
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}
