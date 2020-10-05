[![Build Status](https://travis-ci.com/ssayols/rrvgo.svg?branch=master)](https://travis-ci.com/ssayols/rrvgo)

# RRVGO

Reduce and visualize lists of Gene Ontology terms by identifying redudance based on semantic similarity.

## Introduction to rrvgo

Gene Ontologies (GO) are often used to guide the interpretation of high-throughput
omics experiments, with lists of differentially regulated genes being summarized
into sets of genes with a common functional representation. Due to the hierachical
nature of Gene Ontologies, the resulting lists of enriched sets are usually
redundant and difficult to interpret.

`rrvgo` aims at simplifying the redundance of GO sets by grouping similar terms
in terms of semantic similarity. It also provides some plots to help with
interpreting the summarized terms.

This software is heavily influenced by [REVIGO](http://revigo.irb.hr/). It mimics
a good part of its core functionality, and even some of the outputs are similar.
Without aims to compete, `rrvgo` tries to offer a programatic interface using
available annotation databases and semantic similarity methods implemented in the
Bioconductor project.

## Installation

To install this package, start R and enter:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("rrvgo")
```

## Documentation

To view documentation for the version of this package installed in your system, start R and enter:

```r
browseVignettes("rrvgo")
```

or access the [pkgdown documentation](https://ssayols.github.io/rrvgo/index.html).

## Dependencies

### Imports (mandatory for core functionality)

* [GOSemSim](http://bioconductor.org/packages/GOSemSim/): Semantic similarity computation among GO terms.
* [AnnotationDbi](http://bioconductor.org/packages/AnnotationDbi/): Provides user interface and database connection code for annotation data packages using SQLite data storage.
* [GO.db](http://bioconductor.org/packages/GO.db/): Annotation maps describing the entire Gene Ontology assembled using data from GO.

### Suggests

* pheatmap
* wordcloud
* treemap
* ggplot2
* shiny
* heatmaply
* plotly
* [org.Ag.eg.db](http://bioconductor.org/packages/org.Ag.eg.db/): Genome wide annotation for Anopheles.
* [org.At.tair.db](http://bioconductor.org/packages/org.At.tair.db/): Genome wide annotation for Arabidopsis.
* [org.Bt.eg.db](http://bioconductor.org/packages/org.Bt.eg.db/): Genome wide annotation for Bovine.
* [org.Ce.eg.db](http://bioconductor.org/packages/org.Ce.eg.db/): Genome wide annotation for Worm.
* [org.Cf.eg.db](http://bioconductor.org/packages/org.Cf.eg.db/): Genome wide annotation for Canine.
* [org.Dm.eg.db](http://bioconductor.org/packages/org.Dm.eg.db/): Genome wide annotation for Fly.
* [org.Dr.eg.db](http://bioconductor.org/packages/org.Dr.eg.db/): Genome wide annotation for Zebrafish.
* [org.EcK12.eg.db](http://bioconductor.org/packages/org.EcK12.eg.db/): Genome wide annotation for E coli strain K12.
* [org.EcSakai.eg.db](http://bioconductor.org/packages/org.EcSakai.eg.db/): Genome wide annotation for E coli strain Sakai.
* [org.Gg.eg.db](http://bioconductor.org/packages/org.Gg.eg.db/): Genome wide annotation for Chicken.
* [org.Hs.eg.db](http://bioconductor.org/packages/org.Hs.eg.db/): Genome wide annotation for Human.
* [org.Mm.eg.db](http://bioconductor.org/packages/org.Mm.eg.db/): Genome wide annotation for Mouse.
* [org.Mmu.eg.db](http://bioconductor.org/packages/org.Mmu.eg.db/): Genome wide annotation for Rhesus.
* [org.Pf.plasmo.db](http://bioconductor.org/packages/org.Pf.plasmo.db/): Genome wide annotation for Malaria.
* [org.Pt.eg.db](http://bioconductor.org/packages/org.Pt.eg.db/): Genome wide annotation for Chimp.
* [org.Rn.eg.db](http://bioconductor.org/packages/org.Rn.eg.db/): Genome wide annotation for Rat.
* [org.Sc.sgd.db](http://bioconductor.org/packages/org.Sc.sgd.db/): Genome wide annotation for Yeast.
* [org.Ss.eg.db](http://bioconductor.org/packages/org.Ss.eg.db/): Genome wide annotation for Pig.
* [org.Xl.eg.db](http://bioconductor.org/packages/org.Xl.eg.db/): Genome wide annotation for Xenopus.
