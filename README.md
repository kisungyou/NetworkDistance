
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NetworkDistance

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/NetworkDistance)](https://CRAN.R-project.org/package=NetworkDistance)
[![Travis build
status](https://travis-ci.org/kisungyou/NetworkDistance.svg?branch=master)](https://travis-ci.org/kisungyou/NetworkDistance)
[![](https://cranlogs.r-pkg.org/badges/NetworkDistance)](https://cran.r-project.org/package=NetworkDistance)
<!-- badges: end -->

NetworkDistance package is a collection of **inter-graph** distance
measures. Instead of graph distance that measures the degree of farness
between nodes within a graph, we consider each network as an object and
compute distance between those objects.

## Installation

You can install the released version of NetworkDistance from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("NetworkDistance")
```

or the development version from github:

``` r
## install.packages("devtools")
## library(devtools)
devtools::install_github("kisungyou/NetworkDistance")
```

## Currently Available Methods

We support following methds at this stage and the collection will be
expanded continuously.

| Function        |         Reference          | Description                                               |
|-----------------|:--------------------------:|-----------------------------------------------------------|
| `nd.centrality` |     Roy et al. (2014)      | Distance by Network Centrality Measures                   |
| `nd.csd`        | Ipsen and Mikhailov (2002) | *L*<sub>2</sub> Distance of Continuous Spectral Densities |
| `nd.dsd`        |   Wilson and Zhu (2008)    | Discrete Spectral Distance                                |
| `nd.edd`        |                            | Edge Difference Distance                                  |
| `nd.extremal`   | Jakobson and Rivin (2002)  | Extremal Distance with Top-*k* Eigenvalues                |
| `nd.gdd`        |   Hammond et al. (2013)    | Graph Diffusion Distance                                  |
| `nd.graphon`    |  Mukherjee et al. (2017)   | Graphon Estimates Distance                                |
| `nd.hamming`    |       Hamming (1950)       | Hamming Distance                                          |
| `nd.him`        |    Jurman et al. (2015)    | Hamming-Ipsen-Mikhailov (HIM) Distance                    |
| `nd.moments`    |  Mukherjee et al. (2017)   | Log Moments Distanec                                      |
| `nd.nfd`        |     Bao et al. (2018)      | Network Flow Distance                                     |
| `nd.wsd`        |     Fay et al. (2010)      | Distance with Weighted Spectral Distribution              |
