# BLFRM
Code for implementing Bayesian Linear Factor Regression Model. 

This README describes code that realizes "Bayesian Latent Factor Regression for Functional and Longitudinal Data" by Silvia Montagna, Surya T. Tokdar, Brian Neelon, and David B. Dunson


## Data

The main analysis in the code uses data from processed NHANES accelerometry data. Interested groups can access similar data here: https://github.com/andrew-leroux/rnhanesdata or apply to the repo owner for access. 

This supplement also includes a simulated dataset. The file named "simulated_example_blfrm" generates functional data response using the "dfosr" package; One can access "dfosr" package here .https://rdrr.io/github/drkowal/dfosr/. Code to generate this dataset and others is included in this supplement. 


## Code

The `source` subdirectory contains code used for fitting the BLFRM modelA. There are three files:

* `BLFRM.R` contains the `BLFRM` function, and a description of that function. This is the primary user-facing function. 
* `BLFRM_Utilities.R` contains various helper functions, including data visualization functions, functions for generating basis matrix, and a function to generate the posteriors of multiplicative gamma process (MGP). 

The code has been tested on R version 4.0.2. To install packages required for BLFRM and the code in this supplement, please use:

```{r}
install.packages(c("fda", "coda", "truncdist", "mvtnorm", "dplyr","tidyverse", "dfosr"))
```


## Instructions for use

The `examples` subdirectory contains two reproducible reports with similar structures. These implement the BLFRM model for the datasets described above, plot showing the posterior mean of the factors with the simultaneous and pointwise 95% credible bands, plot showing the posterior mean, simultaneous and pointwise 95% credible bands for a single curve. The simulated data example is smaller, and runs in less than 20 minute, while the example of the processed NHANES accelerometry data runs in roughly one hour. 
