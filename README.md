This directory will contain R scripts for the estimation of the mixing
distribution for k-monotone and completely monotone distributions. We
start by a simulation for 1-monotone distributions, where the mixing
distribution is Poisson(10). The script simulation_1_monotone.R
generates 100 samples of sizes 100, 1000 and 10,000 and computes the
boxplot and draws the estimate of the mixing distribution for the last
sample of size 10,000. The estimate is constructed via the discrete
Grenander estimate, as described in Section 4 of k_monotone.pdf.
The R package Rcpp is used, so Windows users should have Rtools running
and Mac users should have a terminal version of Xcode.
