# StreamVAST
R package for modeling stream networks with VAST

## Description
StreamVAST is a package that aids users designing spatio-temporal models for data collected in stream networks. It is designed to integrate easily with the [tinyVAST package](https://github.com/vast-lib/tinyVAST)), and users will need to be somewhat familiar with tinyVAST to understand and make full use these features. A working knowledge of the [sfnetworks](https://cran.r-project.org/web/packages/sfnetworks/vignettes/sfn01_structure.html) package is also recommended.  
StreamVAST contains a variety of functions to guide users in converting linework from other sources into a valid network with the appropriate characteristics. User have options for selecting a root node, pruning unnecessary branches, dividing the network into prediction frames, and associating various data types with network features. Other features will assist with generating the objects needed for VAST stream network functionality, assessing model fit, and easily making various maps and plots. 

This package is still in the early phases of development, and users should expect frequent updates and changes.
## Installation
It is strongly recommended that users first install the tinyVAST package and its dependencies. 

To install StreamVAST, use this code:
```
remotes::install_github("Jpharris7/StreamVAST")
```
If this doesn't work, consider leaving out the vignettes, as these are under development.
```
remotes::install_github("Jpharris7/StreamVAST", build_vignettes = FALSE)
```
## Shape and Data Prep
This section is demonstrates how to format and clean a set of lines, convert it into a network, root the network, remove unnecessary sections, and associate various types of data with the network. [Preparing a Stream Network](https://jpharris7.github.io/StreamVAST/articles/shape_prep.html)

## Formating and Running VAST
This section shows how to use the network and the associated data to produce a model using VAST functions, and provides some basic advice for choosing settings and other features. [Fitting a VAST Model](https://jpharris7.github.io/StreamVAST/articles/model_fitting.html)
## Plots and outputs
The tutorial for this section is under construction.
