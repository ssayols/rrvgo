#' calcSimilarityScores
#' Calculate the score similarity matrix between terms
#' 
#' @param goterms vector of GO terms
#' @param orgdb one of org.* Bioconductor packages (the package name, or the
#'   package itself)
#' @param semdata object with prepared GO DATA for measuring semantic similarity
#' @param ont ontlogy. One of c("BP", "MF", "CC")
#' @param method distance method. One of the supported methods by GOSemSim:
#'   c("Resnik", "Lin", "Rel", "Jiang", "Wang")
#' @param verbose be explicit while dropping GO terms not found in orgdb
#' @return a square matrix with similarity scores between terms
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example2.txt", package="RVGO"))
#' simMatrix <- calculateSimMatrix(go_analysis$V1, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' @importFrom GOSemSim godata goSim
#' @export
calculateSimMatrix <- function(x, orgdb,
                               semdata=GOSemSim::godata(orgdb, ont=ont),
                               ont=c("BP", "MF", "CC"),
                               method=c("Resnik", "Lin", "Rel", "Jiang", "Wang"),
                               verbose=TRUE) {
 
  # check function args
  ont <- match.arg(ont) 
  method <- match.arg(method) 
  
  # load orgdb object
  if(class(orgdb) != "OrgDb") {
    orgdb <- loadOrgdb(orgdb)
  }
 
  # filter GO terms not in orgdb
  x <- unique(x)
  found <- x %in% names(semdata@IC)
  if(verbose) {
    warning("Removed ", length(x) - sum(found), " terms not found in orgdb for ", ont)
  }
  x <- x[found]
 
  # return the similarity matrix
  matrix(GOSemSim::goSim(x, x, semData=semdata, measure=method),
         ncol=length(x), dimnames=list(x, x))
}

#' reduceSimMatrix
#' Reduce a set of GO terms based on their semantic similarity and scores.
#' 
#' @details
#' Remove terms with a similarity higher than `threshold`. Decide which term
#' remains based on a score. If no score is provided, the select either the
#' broader or the narrower one (`untie` parm).
#' 
#' @param simMatrix a (square) similarity matrix
#' @param scores *named* vector with scores (weights) assigned to each term.
#'   Higher is better. Can be NULL (default, means no scores). Note: if you have
#'   p-values as scores, consider -1*log-transforming them (`-log(p)`)
#' @param threshold similarity threshold
#' @param orgdb one of org.* Bioconductor packages (the package name, or the
#'   package itself)
#' @return a data.frame with all terms and it's "reducer" (NA if the term was not reduced)
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example2.txt", package="RVGO"))
#' scores <- setNames(-log10(go_analysis$V2), go_analysis$V1)
#' simMatrix <- calculateSimMatrix(go_analysis$V1, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' reduced_go_analysis <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' @export
reduceSimMatrix <- function(simMatrix, scores=NULL, threshold, orgdb) {
 
  # check function arguments
  if(!is.null(scores) && !all(rownames(simMatrix) %in% names(scores))) {
    stop("scores vector does not contain all terms in the similarity matrix")
  }
  scores <- scores[match(rownames(simMatrix), names(scores))]
  
  # get category size, and use it as scores if they were not provided
  sizes <- getGoSize(rownames(simMatrix), orgdb)
  if(is.null(scores)) {
    scores <- sizes
  }
  
  # sort matrix based on the score
  orows <- match(rownames(simMatrix), names(scores))
  ocols <- match(colnames(simMatrix), names(scores))  # JIC don't come in order
  simMatrix <- simMatrix[orows, ocols]
  o <- rev(order(scores, sizes, na.last=FALSE))
  simMatrix <- simMatrix[o, o]
  
  # for each term find other terms above similarity threshold
  parent <- character(nrow(simMatrix))
  for(i in (nrow(simMatrix)-1):1) {
    similar <- which(simMatrix[i, (i+1):ncol(simMatrix)] > threshold)
    parent[i+similar] <- rownames(simMatrix)[i]
  }
  
  # return
  data.frame(go=rownames(simMatrix),
             parent=parent,
             scores=scores[match(rownames(simMatrix), names(scores))],
             size=sizes[match(rownames(simMatrix), names(scores))],
             term=getGoTerm(rownames(simMatrix)))
}

#' getGoSize
#' Get GO term size (# of genes)
#' 
#' @param terms GO terms
#' @param orgdb one of org.* Bioconductor packages (the package name, or the
#'   package itself)
#' @importFrom AnnotationDbi select keys
#' @return number of genes associated with each term
getGoSize <- function(terms, orgdb) {
  if(class(orgdb) != "OrgDb") {
    orgdb <- loadOrgdb(orgdb)
  }
  
  # get all GO terms with genes associated
  go <- suppressMessages(
          AnnotationDbi::select(orgdb,
                                keytype="ENTREZID",
                                columns=c("GO", "ONTOLOGY"),
                                keys=AnnotationDbi::keys(orgdb, keytype="ENTREZID")))
  go <- go[!is.na(go$GO), ]
  go <- go[go$GO %in% terms, ]
  
  # count
  counts   <- table(go$GO)
  empty    <- terms[!terms %in% names(counts)]
  nocounts <- setNames(rep(NA, length(empty)), empty)
  
  c(counts, nocounts)
}

#' getGoTerm
#' Get the description of a GO term
#' 
#' @param x GO terms
#' @import GO.db
#' @return the Term slot in GO.db::GOTERM[[x]]
getGoTerm <- function(x) {
  sapply(x, function(x) GO.db::GOTERM[[x]]@Term)
}

#' loadOrgdb
#' Load an orgdb object
#'
#' @param orgdb one of org.* Bioconductor packages
#' @return the loaded orgdb
loadOrgdb <- function(orgdb) {
  if(!requireNamespace(orgdb, quietly = TRUE)) {
    stop("Bioconductor orgdb for", org, "not found. Consider installing it.",
         call. = FALSE)
  }
  eval(parse(text=paste0(orgdb, "::", orgdb)))
}
