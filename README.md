
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
- [Deom on the simulated data](#Deom-on-the-simulated-data)
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
- `exponential_est.stan`: setting up bayesian survival model with the assumption of exponential baseline hazard function;
- `weibull_est.stan`: setting up bayesian survival model with the assumption of weibull baseline hazard function;
- `bSpline_est.stan`:

#### Stan data constructor
- `stan_exponential.R`: setting up the stan_data structure for `exponential_est.stan`;
- `stan_weibull.R`: setting up the stan_data structure for `weibull_est.stan`;
- `stan_bSpline.R`: setting up the stan_data structure for `bSpline_est.stan`.

#### Other supporting functions
- `Functions.R`: Functions being used when constructing the stan_data, model diagnosis and model performance evaluation.

### Report

# Validating the method:
The method was validated from two following perspectivess. 
- `Model diagnosis`: The model diagnosis focused on two levels: (1) MCMC convergence; and (2) model assumption check;

- `Model performance`: The model performance focused on three levels: (1) parameter estimation estimation; (2) survival probability prediction; and (3) variable selection.

# Deom on the simulated data:
- `Sim_data_Analysis.Rmd`:

With the listed code, [a simulated bayesian survival analysis](./Sim_data_Analysis.Rmd) (i.e., a bayesian survival model for time-to -event data under high-dimensional setting was conducted on a simulated dataset) was conducted to validate the research setting. See the model evaluation result on this simulated dataset as following: 

![alt text](./Images/PowerCurve.jpg)



