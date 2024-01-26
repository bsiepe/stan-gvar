## Overview
This repository contains code used to experiment with implementing Bayesian GVAR models in the `tsnet` package. 
It also contains code to reproduce the simulation study comparing the Stan approach to the BGGM approach in the revised version of the preprint "Bayesian Estimation and Comparison of Idiographic Network Models" by Siepe & Heck (2023) (https://psyarxiv.com/uwfjc/). 
Start the .Rproj var-compare.Rproj before running the scripts. 

## Scripts
The file `implementation-comparison-sim.qmd` contains the code for the implementation comparison simulation study. 


## MPlus version
We also implemented a generic MPlus version of Bayesian GVAR models. The wrapper function is included in `scripts\functions.R`. 
An example of its use is given in `GVAR_mplus.qmd`. 

