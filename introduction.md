---
title: "Using RMarkdown"
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do you write a lesson using R Markdown and `{sandpaper}`?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Explain how to use markdown with the new lesson template
- Demonstrate how to include pieces of code, figures, and nested challenge blocks

::::::::::::::::::::::::::::::::::::::::::::::::


## Packages



``` r
# Proteomics


library(limpa)
```

``` output
Loading required package: limma
```

``` r
# scRNA-seq

library(Seurat)
```

``` output
Loading required package: SeuratObject
```

``` output
Loading required package: sp
```

``` output

Attaching package: 'SeuratObject'
```

``` output
The following objects are masked from 'package:base':

    intersect, t
```

``` r
#library(SeuratData)


# RNA-seq pathway

library(clusterProfiler)
```

``` output

```

``` output
clusterProfiler v4.16.0 Learn more at https://yulab-smu.top/contribution-knowledge-mining/

Please cite:

Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
clusterProfiler: an R package for comparing biological themes among
gene clusters. OMICS: A Journal of Integrative Biology. 2012,
16(5):284-287
```

``` output

Attaching package: 'clusterProfiler'
```

``` output
The following object is masked from 'package:stats':

    filter
```

``` r
library(enrichplot)
```

``` output
enrichplot v1.28.4 Learn more at https://yulab-smu.top/contribution-knowledge-mining/

Please cite:

Guangchuang Yu. Using meshes for MeSH term enrichment and semantic
analyses. Bioinformatics. 2018, 34(21):3766-3767,
doi:10.1093/bioinformatics/bty410
```

``` r
library(RegEnrich)
```

``` output
Loading required package: S4Vectors
```

``` output
Loading required package: stats4
```

``` output
Loading required package: BiocGenerics
```

``` output
Loading required package: generics
```

``` output

Attaching package: 'generics'
```

``` output
The following objects are masked from 'package:base':

    as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    setequal, union
```

``` output

Attaching package: 'BiocGenerics'
```

``` output
The following object is masked from 'package:limma':

    plotMA
```

``` output
The following objects are masked from 'package:stats':

    IQR, mad, sd, var, xtabs
```

``` output
The following objects are masked from 'package:base':

    anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    unsplit, which.max, which.min
```

``` output

Attaching package: 'S4Vectors'
```

``` output
The following object is masked from 'package:clusterProfiler':

    rename
```

``` output
The following object is masked from 'package:utils':

    findMatches
```

``` output
The following objects are masked from 'package:base':

    expand.grid, I, unname
```

``` output
Loading required package: dplyr
```

``` output

Attaching package: 'dplyr'
```

``` output
The following objects are masked from 'package:S4Vectors':

    first, intersect, rename, setdiff, setequal, union
```

``` output
The following objects are masked from 'package:BiocGenerics':

    combine, intersect, setdiff, setequal, union
```

``` output
The following object is masked from 'package:generics':

    explain
```

``` output
The following objects are masked from 'package:stats':

    filter, lag
```

``` output
The following objects are masked from 'package:base':

    intersect, setdiff, setequal, union
```

``` output
Loading required package: tibble
```

``` output
Loading required package: BiocSet
```

``` output
Loading required package: SummarizedExperiment
```

``` output
Loading required package: MatrixGenerics
```

``` output
Loading required package: matrixStats
```

``` output

Attaching package: 'matrixStats'
```

``` output
The following object is masked from 'package:dplyr':

    count
```

``` output

Attaching package: 'MatrixGenerics'
```

``` output
The following objects are masked from 'package:matrixStats':

    colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    colWeightedMeans, colWeightedMedians, colWeightedSds,
    colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    rowWeightedSds, rowWeightedVars
```

``` output
Loading required package: GenomicRanges
```

``` output
Loading required package: IRanges
```

``` output

Attaching package: 'IRanges'
```

``` output
The following objects are masked from 'package:dplyr':

    collapse, desc, slice
```

``` output
The following object is masked from 'package:clusterProfiler':

    slice
```

``` output
The following object is masked from 'package:sp':

    %over%
```

``` output
Loading required package: GenomeInfoDb
```

``` output
Loading required package: Biobase
```

``` output
Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.
```

``` output

Attaching package: 'Biobase'
```

``` output
The following object is masked from 'package:MatrixGenerics':

    rowMedians
```

``` output
The following objects are masked from 'package:matrixStats':

    anyMissing, rowMedians
```

``` warning
Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
```

``` output

Attaching package: 'SummarizedExperiment'
```

``` output
The following object is masked from 'package:Seurat':

    Assays
```

``` output
The following object is masked from 'package:SeuratObject':

    Assays
```

``` r
library(STRINGdb)
library(goseq)
```

``` output
Loading required package: BiasedUrn
```

``` output
Loading required package: geneLenDataBase
```

``` output

Attaching package: 'geneLenDataBase'
```

``` output
The following object is masked from 'package:S4Vectors':

    unfactor
```

``` r
library(fgsea)
library(ReactomePA)
```

``` output
ReactomePA v1.52.0 Learn more at https://yulab-smu.top/contribution-knowledge-mining/

Please cite:

Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for
reactome pathway analysis and visualization. Molecular BioSystems.
2016, 12(2):477-479
```

``` r
library(DOSE)
```

``` output
DOSE v4.2.0 Learn more at https://yulab-smu.top/contribution-knowledge-mining/

Please cite:

Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
R/Bioconductor package for Disease Ontology Semantic and Enrichment
analysis. Bioinformatics. 2015, 31(4):608-609
```

``` r
library(pathview)
```

``` output

```

``` output
##############################################################################
Pathview is an open source software package distributed under GNU General
Public License version 3 (GPLv3). Details of GPLv3 is available at
http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
formally cite the original Pathview paper (not just mention it) in publications
or products. For details, do citation("pathview") within R.

The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
license agreement (details at http://www.kegg.jp/kegg/legal.html).
##############################################################################
```

``` r
library(KEGGgraph)
```

``` output

Attaching package: 'KEGGgraph'
```

``` output
The following object is masked from 'package:pathview':

    KEGGEdgeSubtype
```

``` output
The following object is masked from 'package:graphics':

    plot
```

``` output
The following object is masked from 'package:base':

    plot
```

``` r
library(gage)
library(gageData)
library(WGCNA)
```

``` output
Loading required package: dynamicTreeCut
```

``` output
Loading required package: fastcluster
```

``` output

Attaching package: 'fastcluster'
```

``` output
The following object is masked from 'package:stats':

    hclust
```

``` output

Attaching package: 'WGCNA'
```

``` output
The following object is masked from 'package:IRanges':

    cor
```

``` output
The following object is masked from 'package:S4Vectors':

    cor
```

``` output
The following object is masked from 'package:stats':

    cor
```

``` r
library(dplyr)
```





## Introduction

This is a lesson created via The Carpentries Workbench. It is written in
[Pandoc-flavored Markdown](https://pandoc.org/MANUAL.txt) for static files and
[R Markdown][r-markdown] for dynamic files that can render code into output. 
Please refer to the [Introduction to The Carpentries 
Workbench](https://carpentries.github.io/sandpaper-docs/) for full documentation.

What you need to know is that there are three sections required for a valid
Carpentries lesson template:

 1. `questions` are displayed at the beginning of the episode to prime the
    learner for the content.
 2. `objectives` are the learning objectives for an episode displayed with
    the questions.
 3. `keypoints` are displayed at the end of the episode to reinforce the
    objectives.

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: instructor

Inline instructor notes can help inform instructors of timing challenges
associated with the lessons. They appear in the "Instructor View"

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1: Can you do it?

What is the output of this command?

```r
paste("This", "new", "lesson", "looks", "good")
```

:::::::::::::::::::::::: solution 

## Output
 
```output
[1] "This new lesson looks good"
```

:::::::::::::::::::::::::::::::::


## Challenge 2: how do you nest solutions within challenge blocks?

:::::::::::::::::::::::: solution 

You can add a line with at least three colons and a `solution` tag.

:::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::

## Figures

You can also include figures generated from R Markdown:


``` r
pie(
  c(Sky = 78, "Sunny side of pyramid" = 17, "Shady side of pyramid" = 5), 
  init.angle = 315, 
  col = c("deepskyblue", "yellow", "yellow3"), 
  border = FALSE
)
```

<div class="figure" style="text-align: center">
<img src="fig/introduction-rendered-pyramid-1.png" alt="pie chart illusion of a pyramid"  />
<p class="caption">Sun arise each and every morning</p>
</div>

Or you can use standard markdown for static figures with the following syntax:

`![optional caption that appears below the figure](figure url){alt='alt text for
accessibility purposes'}`

![You belong in The Carpentries!](https://raw.githubusercontent.com/carpentries/logo/master/Badge_Carpentries.svg){alt='Blue Carpentries hex person logo with no text.'}

::::::::::::::::::::::::::::::::::::: callout

Callout sections can highlight information.

They are sometimes used to emphasise particularly important points
but are also used in some lessons to present "asides": 
content that is not central to the narrative of the lesson,
e.g. by providing the answer to a commonly-asked question.

::::::::::::::::::::::::::::::::::::::::::::::::


## Math

One of our episodes contains $\LaTeX$ equations when describing how to create
dynamic reports with {knitr}, so we now use mathjax to describe this:

`$\alpha = \dfrac{1}{(1 - \beta)^2}$` becomes: $\alpha = \dfrac{1}{(1 - \beta)^2}$

Cool, right?

::::::::::::::::::::::::::::::::::::: keypoints 

- Use `.md` files for episodes when you want static content
- Use `.Rmd` files for episodes when you need to generate output
- Run `sandpaper::check_lesson()` to identify any issues with your lesson
- Run `sandpaper::build_lesson()` to preview your lesson locally

::::::::::::::::::::::::::::::::::::::::::::::::

[r-markdown]: https://rmarkdown.rstudio.com/
