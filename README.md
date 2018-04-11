# damPassageSims
Simulations for assessing a priori sample sizes for dam passage survival analysis using state-space models

<br>
 
## File descriptions

`singleReleaseBayesian.R` simulation and analysis of fish survival through dams using Cormack-Jolly-Seber mark-recapture survival models. Assumes that all fish are released in the same location.

[JAGS](http://sourceforge.net/projects/mcmc-jags/files/l) must be installed to run the code, and the package `R2jags` must be installed in `R`. The code is modified from [Kéry and Schaub (2010) _Bayesian Population Analysis Using WinBUGS_](http:/www.vogelwarte.ch/de/projekte/publikationen/bpa/).

`post-processing.R` data visualization for simulation results.

## Supporting publications

This simulation was part of a larger analysis:

Zydlewski, J., **D. Stich**, and D. Sigourney. 2017. Hard choices in assessing survival past dams; a comparison of single and paired release strategies. Canadian Journal of Fisheries and Aquatic Sciences 42:178–190. [journal](http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2015-0480#.Ws53wy7wbIU)
