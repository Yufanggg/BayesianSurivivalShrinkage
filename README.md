
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
- [stan data constructor](#eeg-data-analysis)
- [model evaluation](#data-analysis)
- [Project Structure](#project-structure)
- [Results](#Results)

## Requirments
To run this Project, you will need the following:
- R (> 3.6)
<!-- - lmer (install.library("lmer")) 
- lmerTest (install.library(")) --> 

## Installation

## Experimental Design

### Power Analysis
With a given experimental materials, a [Power Analysis](./DOE.Rmd) was conducted to validate the number of participants. See the an example result for 2-by-2 experimental design within subject and target word (which are fully crossed) as following: 
![alt text](./Images/PowerCurve.jpg)

## Data Analysis
### Behavioral Data Analysis
The behavioral data analysis can be conducted:
- `Exp02_BehaviouralAnalysis.Rmd`:

### EEG Data Analysis
EEG data analysis includes: data preprocessing, feature extraction and statistcial modelling.

##### EEG Data Preprocessing \& Feature Extraction can be conducted with the code in the given steps: 
- `EEG_Step01_PreproEEG_Batch.m`:
- `EEG_Step02_PreproEEG.m`:
- `EEG_Step03_PreproRemovthreshold.m`:
- `EEG_Step_04_ExploratoryERPPermutation.m`:

##### EEG statistcial modelling can be conducted with the code in the given steps: 
- `Exp02_ExploratoryERP.Rmd`:
- `Exp02_ERPPlannedAnalysis.Rmd`:
