
# Cognitive2Computation
[![LinkedIn](https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555)](https://www.linkedin.com/in/yufang-w-1295881b5/) [![GitHub](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white&colorB=555)](https://github.com/Yufanggg) <img alt="GitHub" src="https://img.shields.io/github/license/bopith/UnicornCompanies?style=for-the-badge"> 

## Overview
This repository contains the bayesian survival model for high-dimensional data. 

This project is updated from codes used in my Master work and aims to address the following research related questions:
1. The reformulated \textit{Bayint} model for time-to-event data;
2. Its potential pros and cons compared to tree-based machine learning methods;
3. An effective approach to estimate interaction effects for time-to-event data with a moderate n to p ratio.


## Table of Contents

- [stan model](#stan-model)
- [stan data constructor](#stan-data-constructor)
- [model evaluation](#data-analysis)
- [Project Structure](#project-structure)
- [Results](#Results)

## Requirments
To run this Project, you will need the following:
- R (> 3.6)
<!-- - lmer (install.library("lmer")) 
- lmerTest (install.library(")) --> 

## Installation

## Project structure
### Code
#### Stan model:
- `exponential_est.stan`:
- `weibull_est.stan`:
- `bSpline_est.stan`:

#### Stan data constructor
- `stan_exponential.R`:
- `stan_weibull.R`:
- `stan_bSpline.R`:

#### Other supporting functions
- `Functions.R`:

# Deom on the simulated data:
- `Sim_data_Analysis.Rmd`:

With the listed code, [a simulated bayesian survival analysis](./Sim_data_Analysis.Rmd) (i.e., a bayesian survival model for time-to -event data under high-dimensional setting was conducted on a simulated dataset) was conducted to validate the research setting. See the model evaluation result on this simulated dataset as following: 

![alt text](./Images/PowerCurve.jpg)



