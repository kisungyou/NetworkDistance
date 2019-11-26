
<!-- README.md is generated from README.Rmd. Please edit that file -->
NetworkDistance
===============

<!-- badges: start -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/NetworkDistance?color=green)](https://cran.r-project.org/package=NetworkDistance) [![Travis build status](https://travis-ci.org/kyoustat/NetworkDistance.svg?branch=master)](https://travis-ci.org/kyoustat/NetworkDistance) [![](https://cranlogs.r-pkg.org/badges/NetworkDistance)](https://cran.r-project.org/package=NetworkDistance) <!-- badges: end -->

NetworkDistance package is a collection of **inter-graph** distance measures. Instead of graph distance that measures the degree of farness between nodes within a graph, we consider each network as an object and compute distance between those objects.

Installation
------------

You can install the released version of NetworkDistance from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("NetworkDistance")
```

or the development version from github:

``` r
## install.packages("devtools")
## library(devtools)
devtools::install_github("kyoustat/NetworkDistance")
```

Currently Available Methods
---------------------------

We support following methds at this stage and the collection will be expanded continuously.

<table>
<colgroup>
<col width="16%" />
<col width="30%" />
<col width="52%" />
</colgroup>
<thead>
<tr class="header">
<th>Function</th>
<th align="center">Reference</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>nd.centrality</code></td>
<td align="center">Roy et al. (2014)</td>
<td>Distance by Network Centrality Measures</td>
</tr>
<tr class="even">
<td><code>nd.csd</code></td>
<td align="center">Ipsen and Mikhailov (2002)</td>
<td><span class="math inline"><em>L</em><sub>2</sub></span> Distance of Continuous Spectral Densities</td>
</tr>
<tr class="odd">
<td><code>nd.dsd</code></td>
<td align="center">Wilson and Zhu (2008)</td>
<td>Discrete Spectral Distance</td>
</tr>
<tr class="even">
<td><code>nd.edd</code></td>
<td align="center"></td>
<td>Edge Difference Distance</td>
</tr>
<tr class="odd">
<td><code>nd.extremal</code></td>
<td align="center">Jakobson and Rivin (2002)</td>
<td>Extremal Distance with Top-<span class="math inline"><em>k</em></span> Eigenvalues</td>
</tr>
<tr class="even">
<td><code>nd.gdd</code></td>
<td align="center">Hammond et al. (2013)</td>
<td>Graph Diffusion Distance</td>
</tr>
<tr class="odd">
<td><code>nd.hamming</code></td>
<td align="center">Hamming (1950)</td>
<td>Hamming Distance</td>
</tr>
<tr class="even">
<td><code>nd.him</code></td>
<td align="center">Jurman et al. (2015)</td>
<td>Hamming-Ipsen-Mikhailov (HIM) Distance</td>
</tr>
<tr class="odd">
<td><code>nd.nfd</code></td>
<td align="center">Bao et al. (2018)</td>
<td>Network Flow Distance</td>
</tr>
<tr class="even">
<td><code>nd.wsd</code></td>
<td align="center">Fay et al. (2010)</td>
<td>Distance with Weighted Spectral Distribution</td>
</tr>
</tbody>
</table>
