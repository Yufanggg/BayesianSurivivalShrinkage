
# Cognitive2Computation
[![LinkedIn](https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555)](https://www.linkedin.com/in/yufang-w-1295881b5/) [![GitHub](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white&colorB=555)](https://github.com/Yufanggg) <img alt="GitHub" src="https://img.shields.io/github/license/bopith/UnicornCompanies?style=for-the-badge"> 

## Overview
This repository contains the bayesian survival model for high-dimensional data. 

This project is updated from codes used in my Master work and aims to address the following research related questions:
1. The reformulated \textit{Bayint} model for time-to-event data;
2. Its potential pros and cons compared to tree-based machine learning methods;
3. An effective approach to estimate interaction effects for time-to-event data with a moderate n to p ratio.


## Table of Contents

- [Requirments](#Requirments)
- [Installation](#Installation)
- [Project structure](#Project-structure)
- [Resuls on the simulated data](#Results-on-the-simulated-data)
- [Results on the real-world data](#Results-on-the-real-world-data)

## Requirments
To run this Project, you will need the following:
- R (> 3.6)
<!-- - lmer (install.library("lmer")) 
- lmerTest (install.library(")) --> 

## Installation

## Project structure 
### Code
#### Stan model:
- `.\stan_files\exponential_est.stan`: setting up bayesian survival model with the assumption of exponential baseline hazard function;
- `.\stan_files\weibull_est.stan`: setting up bayesian survival model with the assumption of weibull baseline hazard function;
- `.\stan_files\bSpline_est.stan`:

#### Stan data constructor
- `.\function_file\exponential_stan_constructor.R`: setting up the stan_data structure for `exponential_est.stan`;
- `.\function_file\weibull_stan_constructor.R`: setting up the stan_data structure for `weibull_est.stan`;
- `.\function_file\bSpline_stan_constructor.R`: setting up the stan_data structure for `bSpline_est.stan`.

#### Other supporting functions
- `.\function_file\Functions.R`: Functions being used when constructing the stan_data, model diagnosis and model performance evaluation.

### Report
The report of this project is now underdrafting.

## Validating the method:
The method was validated from two following perspectivess. 
1. `Model diagnosis`: The model diagnosis focused on two levels: (1) MCMC convergence; and (2) model assumption check;

2. `Model performance`: The model performance focused on three levels: (1) parameter estimation estimation; (2) survival probability prediction; and (3) variable selection.

## Results on the simulated data:
With the listed code, [a simulated bayesian survival analysis](./Sim_data_Analysis.Rmd) (i.e., a bayesian survival model for time-to -event data under high-dimensional setting was conducted on a simulated dataset) was conducted to validate the research setting. The simulated data can be obtained from the file of [Simulated_Data](./Data/imputed_NOTR_DGF.rds), in which a weibull baseline hazard function were assumed. With this simulated data, the following three models have been built: 
1. **the exponential model**: ![alt text](./Image/exponential_performance%20_simulated.png)
2. **the weibull model**:  ![alt text](./Image/weibull_performance%20_simulated.png)
3. **the bSpline models**: ![alt text](./Image/bSpline_performance%20_simulated.png)



## Results on the real-world data:
The propsed method was also [validated on real-world data](./Real_data_Analysis.Rmd). See the model evaluation results as following:

