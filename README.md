# SAR-SAS
Before SAS/ETS 14.2, there is no an intended process for estimating spatial regressions. I proposed a series of three macros to do that. The macros was presented during the 22nd day-conference for SAS users at Quebec city. The results obtained are also confirmed by comparing them with those provided by R packages.

You can get my the presentation [[in french](https://ahoundetoungan.github.io/files/SAS%20SPAT%20REG/SAS_QC_22nd_day_beamer.pdf)] or on [SAS users website](http://wp.clubsasquebec.ca/).
I provided three macros: the first for Moran and Geary tests, the second for spatial regressions estimation and the last for specification tests. Get the macros files below:

[Moran and Geary tests](https://github.com/ahoundetoungan/SpatRegSAS/blob/master/macros/MoranGeary.sas) – [Spatial regressions estimation](https://github.com/ahoundetoungan/SpatRegSAS/blob/master/macros/Spatialregressions.sas) – [Specification tests](https://github.com/ahoundetoungan/SpatRegSAS/blob/master/macros/Specificationtests.sas)

Here is an example with simulated data on africa’s countries. The interest variable is agricultural production, modeled through the temperature, the precipitations and the tax on production.  This is the dataset used and the weights matrix.

Get the example code and data [here](https://github.com/ahoundetoungan/SpatRegSAS/tree/master/example).
