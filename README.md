# StreamVAST
R package for modeling stream networks with VAST

## Description
StreamVAST is a package that aids users designing spatio-temporal models for data collected in stream networks. It is designed to integrate easily with the [VAST package](https://github.com/James-Thorson-NOAA/VAST), and users will need to be familiar with VAST to understand and make full use these features. A working knowledge of the [sfnetworks](https://cran.r-project.org/web/packages/sfnetworks/vignettes/sfn01_structure.html) package is also recommended.  
StreamVAST contains a variety of functions to guide users in converting linework from other sources into a valid network with the appropriate characteristics. User have options for selecting a root node, pruning unnecessary branches, dividing the network into prediction frames, and associating various data types with network features. Other features will assist with generating the objects needed for VAST stream network functionality, assessing model fit, and easily making various maps and plots. 

This package is still in the early phases of development, and users should expect frequent updates and changes.
## Installation
It is strongly recommended that users first install the VAST package and its dependencies. 
```
remotes::install_github("Jpharris7/StreamVAST")
```
If this doesn't work, consider leaving out the vignettes, as these are under development.
```
remotes::install_github("Jpharris7/StreamVAST", build_vignettes = FALSE)
```
## Shape and Data Prep
This section is nearly complete and will be available shortly.
## Formating and Running VAST
The tutorial for this section is under construction.
## Plots and outputs
The tutorial for this section is under construction.
