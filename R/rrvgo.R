#' calculateSimMatrix
#' Calculate the score similarity matrix between terms
#' 
#' @param x vector of GO terms
#' @param orgdb one of org.* Bioconductor packages (the package name, or the
#'   package itself)
#' @param keytype keytype passed to AnnotationDbi::keys to retrieve GO terms 
#'   associated to gene ids in your orgdb
#' @param semdata object with prepared GO DATA for measuring semantic similarity
#' @param ont ontlogy. One of c("BP", "MF", "CC")
#' @param method distance method. One of the supported methods by GOSemSim:
#'   c("Resnik", "Lin", "Rel", "Jiang", "Wang")
#' @return a square matrix with similarity scores between terms
#' @details 
#' All similarity measures available are those implemented in the
#' [GOSemSim package](https://www.bioconductor.org/packages/release/bioc/html/GOSemSim.html),
#' namely the Resnik, Lin, Relevance, Jiang and Wang methods. See the
#' [Semantic Similarity Measurement Based on GO](https://www.bioconductor.org/packages/release/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html#semantic-similarity-measurement-based-on-go)
#' section from the GOSeSim documentation for more details.
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
#' simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' @import GOSemSim
#' @importFrom methods is
#' @export
calculateSimMatrix <- function(x,
                               orgdb,
                               keytype="ENTREZID",
                               semdata=GOSemSim::godata(orgdb, ont=ont, keytype=keytype),
                               ont=c("BP", "MF", "CC"),
                               method=c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
 
  # check function args
  ont <- match.arg(ont) 
  method <- match.arg(method) 
  
  # load orgdb object
  if(all(is(orgdb) != "OrgDb")) {
    orgdb <- loadOrgdb(orgdb)
  }
 
  # filter GO terms not in orgdb
  x <- unique(x)
  found <- x %in% names(semdata@IC)
  hasAncestor <- !is.na(sapply(x, function(x) tryCatch(GOSemSim:::getAncestors(ont)[x], error=function(e) NA)))
  if(all(!found)) {
    warning("No terms were found in orgdb for ", ont,
            "\nCheck that the organism and ontology match the ones provided by orgdb")
    return(NA)
  } else if(!all(found)) {
    warning("Removed ", length(x) - sum(found), " terms that were not found in orgdb for ", ont)
  }
  x <- x[found & hasAncestor]
  
 
  # return the similarity matrix
  m <- matrix(GOSemSim::goSim(x, x, semData=semdata, measure=method),
              ncol=length(x), dimnames=list(x, x))
  
  # removing terms which the similarity couldn't be calculated
  out <- apply(m, 2, function(x) all(is.na(x)))
  m[!out, !out]
}

#' reduceSimMatrix
#' Reduce a set of GO terms based on their semantic similarity and scores.
#' 
#' @details
#' Currently, rrvgo uses the similarity between pairs of terms to compute a 
#' distance matrix, defined as (1-simMatrix). The terms are then hierarchically 
#' clustered using complete linkage, and the tree is cut at the desired threshold, 
#' picking the term with the highest score as the representative of each group.  
#' 
#' Therefore, higher thresholds lead to fewer groups, and the threshold should be 
#' read as the expected similarity of terms within a group (though this is not 
#' entirely correct, and you'll see similarities below this threshold being put in 
#' the same group).
#' 
#' @param simMatrix a (square) similarity matrix
#' @param scores *named* vector with scores (weights) assigned to each term.
#'   Higher is better. Can be NULL (default, means no scores. In this case, a default score
#'   based on set size is assigned, thus favoring larger sets). Note: if you have
#'   p-values as scores, consider -1*log-transforming them (`-log(p)`)
#' @param threshold similarity threshold (0-1). Some guidance:
#'  Large (allowed similarity=0.9), Medium (0.7), Small (0.5), Tiny (0.4)
#'  Defaults to Medium (0.7)
#' @param orgdb one of org.* Bioconductor packages (the package name, or the
#'   orgdb object itself)
#' @param keytype keytype passed to AnnotationDbi::keys to retrieve GO terms 
#'   associated to gene ids in your orgdb
#' @return a data.frame with all terms and it's "reducer" (NA if the term was not reduced)
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
#' simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
#' reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' @importFrom stats cutree hclust
#' @export
reduceSimMatrix <- function(simMatrix, scores=NULL, threshold=0.7, orgdb, keytype="ENTREZID") {
 
  # check function arguments
  if(!is.null(scores) && !all(rownames(simMatrix) %in% names(scores))) {
    stop("Scores vector does not contain all terms in the similarity matrix")
  }

  # get category size, and use it as scores if they were not provided
  sizes <- getGoSize(rownames(simMatrix), orgdb, keytype)
  if(is.null(scores)) {
    message("No scores provided. Falling back to term's size")
    scores <- sizes
  }
  
  scores <- scores[match(rownames(simMatrix), names(scores))]
  
  # reorder the similarity matrix as in the scores, just in case they don't come in the same order
  orows <- match(rownames(simMatrix), names(scores))
  ocols <- match(colnames(simMatrix), names(scores))
  simMatrix <- simMatrix[orows, ocols]
  
  # sort matrix based on the score
  o <- rev(order(scores, sizes, na.last=FALSE))
  simMatrix <- simMatrix[o, o]
  
  # cluster terms and cut the tree at the desired threshold.
  # Then find the term with the highest score as the representative of each cluster
  cluster <- cutree(hclust(as.dist(1-simMatrix)), h=threshold)
  clusterRep <- tapply(rownames(simMatrix), cluster, function(x) x[which.max(scores[x])])
  
  # return
  data.frame(go=rownames(simMatrix),
             cluster=cluster,
             parent=clusterRep[cluster],
             parentSimScore=unlist(Map(seq_len(nrow(simMatrix)), clusterRep[cluster], f=function(i, j) simMatrix[i, j])),
             score=scores[match(rownames(simMatrix), names(scores))],
             size=sizes[match(rownames(simMatrix), names(sizes))],
             term=getGoTerm(rownames(simMatrix)),
             parentTerm=getGoTerm(clusterRep[cluster]))
}

#' getGoSize
#' Get GO term size (# of genes)
#' 
#' @param terms GO terms
#' @param orgdb one of org.* Bioconductor packages (the package name, or the
#'   package itself)
#' @param keytype keytype passed to AnnotationDbi::keys to retrieve GO terms 
#'   associated to gene ids in your orgdb
#' @importFrom AnnotationDbi select keys
#' @importFrom stats setNames
#' @importFrom methods is
#' @return number of genes associated with each term
getGoSize <- function(terms, orgdb, keytype) {
  if(all(is(orgdb) != "OrgDb")) {
    orgdb <- loadOrgdb(orgdb)
  }
  
  # get all GO terms with genes associated
  go <- suppressMessages(
          AnnotationDbi::select(orgdb,
                                keytype=keytype,
                                columns=c("GO", "ONTOLOGY"),
                                keys=AnnotationDbi::keys(orgdb, keytype=keytype)))
  go <- go[!is.na(go$GO), ]
  go <- go[go$GO %in% terms, ]
  
  # count
  counts   <- table(go$GO)
  go <- go[go$GO %in% terms, ]
  empty    <- terms[!(terms %in% names(counts))]
  nocounts <- setNames(rep(0, length(empty)), empty)
  
  c(counts, nocounts)
}

#' getGoTerm
#' Get the description of a GO term
#' 
#' @param x GO terms
#' @import GO.db
#' @return the Term slot in GO.db::GOTERM[[x]]
getGoTerm <- function(x) {
  sapply(x, function(x) tryCatch(GO.db::GOTERM[[x]]@Term, error=function(e) NA))
}

#' loadOrgdb
#' Load an orgdb object
#'
#' @param orgdb one of org.* Bioconductor packages
#' @return the loaded orgdb
loadOrgdb <- function(orgdb) {
  if(!requireNamespace(orgdb, quietly = TRUE)) {
    stop("Bioconductor orgdb for ", orgdb, " not found. Consider installing it.",
         call. = FALSE)
  }
  eval(parse(text=paste0(orgdb, "::", orgdb)))
}
