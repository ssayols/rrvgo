#' scatterPlot
#' Plot GO terms as scattered points.
#' 
#' @param simMatrix a (square) similarity matrix.
#' @param reducedTerms a data.frame with the reduced terms from reduceSimMatrix()
#' @param label add labels with the most representative term of the group.
#' @param pal use custom palette. Defaults to ggplot2's default categorical pal.
#' @details  Distances between points represent the similarity between terms.
#' Axes are the first 2 components of applying a PCA to the similarity matrix. 
scatterPlot <- function(simMatrix, reducedTerms, label=FALSE, pal=NULL) {
  
}

#' treemapPlot
#' Plot GO terms as a treemap.
#' 
#' @param simMatrix a (square) similarity matrix.
#' @param reducedTerms a data.frame with the reduced terms from reduceSimMatrix()
#' @param label add labels with the most representative term of the group.
#' @param pal use custom palette. Defaults to ggplot2's default categorical pal.
#' @details  Distances between points represent the similarity between terms.
#' Axes are the first 2 components of applying a PCA to the similarity matrix.
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example2.txt", package="RVGO"))
#' simMatrix <- calculateSimMatrix(go_analysis$V1, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' reduced_go_analysis <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' treemapPlot(simMatrix, reduced_go_analysis)
treemapPlot <- function(simMatrix, reducedTerms, label=FALSE, pal=NULL) {

}

#' wordlcoudPlot
#' Plot GO reduced terms as a wordcloud.
#' 
#' @param reducedTerms a data.frame with the reduced terms from reduceSimMatrix()
#' @param pal use custom palette. Defaults to "black"
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example2.txt", package="RVGO"))
#' scores <- setNames(-log10(go_analysis$V2), go_analysis$V1)
#' simMatrix <- calculateSimMatrix(go_analysis$V1, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' reduced_go_analysis <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' wordcloudPlot(reduced_go_analysis)
wordcloudPlot <- function(reducedTerms, pal="black") {
  if(!all(sapply(c("wordcloud", "tm", "slam"), requireNamespace, quietly=TRUE))) {
    stop("Package wordcloud and/or its dependencies (tm, slam) not available. ",
         "Consider installing it before using this function.", call.=FALSE)
  }
    
  x <- tm::Corpus(tm::VectorSource(reducedTerms$term))
  x <- tm::tm_map(x, tm::PlainTextDocument)
  x <- tm::tm_map(x, tm::stripWhitespace)
  x <- tm::tm_map(x, tm::removePunctuation)
  x <- tm::tm_map(x, function(x) tm::removeWords(x, tm::stopwords()))
  tdm <- tm::TermDocumentMatrix(x)
  m <- as.matrix(tdm)
  v <- sort(rowSums(m), decreasing=TRUE)
  d <- data.frame(word=names(v), freq=v)
  wordcloud::wordcloud(d$word, d$freq, min.freq=2, colors=pal)
}
