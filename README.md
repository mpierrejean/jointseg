# jointseg

This package implements functions to quickly segment multivariate signals into piecewise constant profiles and a framework to generate realistic copy-number profiles. A typical application is the joint segmentation of total DNA copy numbers and allelic ratios obtained from Single Nucleotide Polymorphism (SNP) microarrays in cancer studies.

# Installation

The package can be installed from github:

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("mpierrejean/jointseg")
```

# Examples and descriptions
A detailled example can be found in [PSSeg.pdf](https://github.com/mpierrejean/jointseg/blob/master/vignettes/PSSeg.pdf) and [dataGeneration.pdf](https://github.com/mpierrejean/jointseg/blob/master/vignettes/dataGeneration.pdf)


## Software status

| Resource:     | GitHub        | Travis CI      | Appveyor         |
| ------------- | ------------------- | -------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & OS X_ | _Windows_        |
| R CMD check  | | <a href="https://travis-ci.org/mpierrejean/jointseg"><img src="https://travis-ci.org/mpierrejean/jointseg.svg" alt="Build status"></a> | <a href="https://ci.appveyor.com/project/mpierrejean/jointseg"><img src="https://ci.appveyor.com/api/projects/status/github/mpierrejean/jointseg?svg=true" alt="Build status"></a> |
| Test coverage | | <a href="https://codecov.io/gh/mpierrejean/jointseg"><img src="https://codecov.io/gh/mpierrejean/jointseg/branch/master/graph/badge.svg" alt="Coverage Status"/></a> | |
