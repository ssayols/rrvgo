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
#' Group terms which are at least within a similarity below `threshold`. Decide 
#' which term remains based on a score. If no score is provided, then decide based 
#' on the "uniqueness" or the term "size".
#'
#' Currently, rrvgo uses the similarity between pairs of terms to compute a 
#' distance matrix, defined as (1-simMatrix). The terms are then hierarchically 
#' clustered using complete linkage, and the tree is cut at the desired threshold, 
#' picking the term with the highest score as the representative of each group.  
#' 
#' Therefore, higher thresholds lead to fewer groups, and the threshold should be 
#' read as the minimum similarity between group representatives.
#' 
#' @param simMatrix a (square) similarity matrix
#' @param scores one of c("uniqueness", "size"), or a *named* vector with scores
#'   provided for each term, where higher values favor choosing the term as the
#'   cluster representative. The default "uniqueness" uses a score reflecting how
#'   unique the term is. Note: if you like to use p-values as scores, consider
#'   -1*log-transforming them (`-log(p)`)
#' @param threshold similarity threshold (0-1). Some guidance:
#'  Large (allowed similarity=0.9), Medium (0.7), Small (0.5), Tiny (0.4)
#'  Defaults to Medium (0.7)
#' @param orgdb one of org.* Bioconductor packages (the package name, or the
#'   orgdb object itself)
#' @param keytype keytype passed to AnnotationDbi::keys to retrieve GO terms 
#'   associated to gene ids in your orgdb
#' @param children when retrieving GO term size, include genes in children terms.
#'   (based on relationships in the GO DAG hierarchy). Defaults to TRUE
#' @return a data.frame identifying the different clusters of terms, the parent
#' term representing the cluster, and some metrics of importance describing how
#' unique and dispensable a term is.
#' @examples
#' go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
#' simMatrix <- calculateSimMatrix(go_analysis$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
#' scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
#' reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
#' @importFrom stats cutree hclust
#' @export
reduceSimMatrix <- function(simMatrix, scores=c("uniqueness", "size"),
                            threshold=0.7, orgdb, keytype="ENTREZID", children=TRUE) {
 
  # check function arguments
  if(is(scores, "character")) {
    scores <- match.arg(scores)
  } else {
    stopifnot("Scores vector does not contain all terms in the similarity matrix" = all(rownames(simMatrix) %in% names(scores)))
  }

  # cluster terms and cut the tree at the desired threshold.
  cluster <- cutree(hclust(as.dist(1 - simMatrix)), h=threshold)

  # get category size and term uniqueness, and use it as scores if they were not provided
  sizes    <- getGoSize(rownames(simMatrix), orgdb, keytype, children)
  termUniq <- getTermUniq(simMatrix, cluster)
  if(is(scores, "character")) {
    scores <- switch(scores,
                     uniqueness={ message("No scores provided. Falling back to term's uniqueness"); termUniq },
                     size      ={ message("No scores provided. Falling back to term's GO size"); sizes })
  }
  scores <- scores[match(rownames(simMatrix), names(scores))]   # shortlist + order scores for terms present in the simMatrix
  
  # sort everything based on the score
  o <- rev(order(scores, termUniq, na.last=FALSE))
  scores    <- scores[o]
  simMatrix <- simMatrix[o, o]
  
  # finally, find the term with the highest score as the representative of each cluster
  cluster    <- cluster[match(rownames(simMatrix), names(cluster))]   # reorder the cluster vector to match row order in the simMatrix
  clusterRep <- tapply(rownames(simMatrix), cluster, function(x) x[which.max(scores[x])])
  
  # return
  data.frame(go                =rownames(simMatrix),
             cluster           =cluster,
             parent            =clusterRep[cluster],
             score             =scores[match(rownames(simMatrix), names(scores))],
             size              =sizes[match(rownames(simMatrix), names(sizes))],
             term              =getGoTerm(rownames(simMatrix)),
             parentTerm        =getGoTerm(clusterRep[cluster]),
             termUniqueness    =getTermUniq(simMatrix),
             termUniquenessWithinCluster=getTermUniq(simMatrix, cluster),
             termDispensability=getTermDisp(simMatrix, cluster, clusterRep))
}

#' getTermUniq
#' Calculate the term uniqueness score, defined as 1 minus the average semantic
#' similarity of a term to all other terms.
#' 
#' @param simMatrix a (square) similarity matrix
#' @param cluster vector with the cluster each entry in the simMatrix belongs to.
#'    If NULL, a 
#' @return a vector of term uniqueness scores
getTermUniq <- function(simMatrix, cluster=NULL) {
  diag(simMatrix) <- 0    # parent terms are completely unique (considered not similar to themselves)
  if(is.null(cluster)) {  # global uniqueness, comparing to all terms regardless of the cluster they belong to
    1 - apply(simMatrix, 1, mean, na.rm=TRUE)
  } else {                # local uniqueness within the cluster
    x <- tapply(names(cluster), cluster, function(x) getTermUniq(simMatrix[x, x, drop=FALSE]), simplify=FALSE)
    names(x) <- NULL
    x <- unlist(x)
    x[rownames(simMatrix)]
  }
}

#' getTermDisp
#' Calculate the term dispensability score, defined as the semantic similarity
#' threshold a term was assigned to a cluster (namely, the similarity of a term
#' to the cluster representative term).
#' 
#' @param simMatrix a (square) similarity matrix
#' @param cluster the cluster assignment for each term
#' @param clusterRep the cluster representative term
#' @return a vector of term dispensability scores
getTermDisp <- function(simMatrix, cluster, clusterRep) {
  unlist(Map(rownames(simMatrix), cluster, f=function(term, cluster) {
    if(term != clusterRep[cluster]) {
      simMatrix[term, clusterRep[cluster]]
    } else {
      0       # parent terms are not dispensable
    }
  }))
}


#' getGoSize
#' Get GO term size (# of genes)
#' 
#' @param terms GO terms
#' @param orgdb one of org.* Bioconductor packages (the package name, or the
#'   package itself)
#' @param keytype keytype passed to AnnotationDbi::keys to retrieve GO terms 
#'   associated to gene ids in your orgdb
#' @param children include genes in children terms (based on relationships in
#'  the GO DAG hierarchy)
#' @importFrom AnnotationDbi select keys
#' @importFrom stats setNames
#' @importFrom methods is
#' @return number of genes associated with each term
getGoSize <- function(terms, orgdb, keytype, children) {
  if(all(is(orgdb) != "OrgDb")) {
    orgdb <- loadOrgdb(orgdb)
  }
  
  # retrieve genes per term and count unique genes within each term
  go <- AnnotationDbi::select(orgdb, keytype=if(children) "GOALL" else "GO", keys=terms, columns=keytype)
  
  # count
  counts   <- tapply(go[[keytype]], go[[if(children) "GOALL" else "GO"]], function(x) length(unique(x)))
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
