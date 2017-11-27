## jointseg

This package implements functions to quickly segment multivariate signals into piecewise-constant profiles, as well as a framework to generate realistic copy-number profiles. A typical application is the joint segmentation of total DNA copy numbers and allelic ratios obtained from Single Nucleotide Polymorphism (SNP) microarrays in cancer studies.

## Installation

You can install jointseg from github with:

    # install.packages("devtools")
    devtools::install_github("mpierrejean/jointseg")

## Usage

The main high-level joint segmentation functions are:
* `jointSeg` for arbitrary signals, see `?jointSeg`.
* `PSSeg` for bivariate copy-number signals, see `?PSSeg` and `vignette("PSSeg")`.

We also refer to  `vignette("dataGeneration")` for a description of the generation of synthetic DNA copy-number profiles using data from the `acnr` package.

## References

Pierre-Jean, M, Rigaill, G. J. and Neuvial, P. (2015). "Performance Evaluation of DNA Copy Number Segmentation Methods." *Briefings in Bioinformatics*, no. 4: 600–615.


## Software status

| Resource:     | GitHub        | Travis CI      | Appveyor         |
| ------------- | ------------------- | -------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & OS X_ | _Windows_        |
| R CMD check  | | <a href="https://travis-ci.org/mpierrejean/jointseg"><img src="https://travis-ci.org/mpierrejean/jointseg.svg" alt="Build status"></a> | <a href="https://ci.appveyor.com/project/mpierrejean/jointseg"><img src="https://ci.appveyor.com/api/projects/status/github/mpierrejean/jointseg?svg=true" alt="Build status"></a> |
| Test coverage | | <a href="https://codecov.io/gh/mpierrejean/jointseg"><img src="https://codecov.io/gh/mpierrejean/jointseg/branch/master/graph/badge.svg" alt="Coverage Status"/></a> | |
