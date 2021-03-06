# BLFRM
An R package for implementing Bayesian Latent Factor Regression Model (BLFRM). 

One can learn more about BLFRM in the paper "Bayesian Latent Factor Regression for Functional and Longitudinal Data" by Silvia Montagna, Surya T. Tokdar, Brian Neelon, and David B. Dunson.


## Data

The main analysis in the code uses processed COVID19 fatality count data (last updated 12/02/2020). The data can be found in the "data" directory. Details for processing the data can be found in the script "source/process_covid_data.R". 
This supplement also includes a simulated dataset. The file named "simulated_example_blfrm" generates functional data response using the "dfosr" package; One can learn more about the "dfosr" package at https://rdrr.io/github/drkowal/dfosr/. Code to generate this dataset and visualize posterior distribution is included in this package. 


## Code

The `source` subdirectory contains code used for fitting the BLFRM. There are three files:

* `BLFRM.R` contains the `BLFRM` function and a description of that function. This is the primary user-facing function. 
* `BLFRM_Utilities.R` contains various helper functions, including data visualization functions, functions for generating basis matrix, and a function to generate the posteriors of multiplicative gamma process (MGP). 
* `process_covid.R` contains code for processing COVID19, Population and Mobility data, and organizing them for analysis.
The code has been tested on R version 4.0.2. To install packages required for BLFRM and the code in this supplement, please use:

```{r}
install.packages(c("fda", "coda", "truncdist", "mvtnorm", "dplyr","tidyverse")
install_github('drkowal/dfosr')
```


## Instructions for use

The `examples` subdirectory contains two reproducible reports with similar structures. These implement the BLFRM for the datasets described above, create a plot showing the posterior mean, simultaneous and pointwise 95% credible bands for a single curve, and the results for effective sample size(Eff). The simulated data example is smaller and runs in less than 20 minutes, while the example of the processed COVID accelerometry data runs longer than one hour. 
