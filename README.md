# DICoP Anomaly Detection Presentation
Annual EMS and Trauma Registry Anomaly Detection Using Julia and Quarto RevealJS

## Overview
This project provides the analytic workflow, statistical methodology, and presentation materials for detecting anomalous reporting patterns in the Iowa EMS and Trauma Registry. The work integrates Julia for data processing, Tidier for data manipulation, and RevealJS via Quarto for presentation delivery.

## Data Description
The analysis uses annual record submission counts for each EMS agency and verified trauma center across Iowa from 2020 through 2026. The data represent stable, operational reporting processes with minimal short‑term volatility. Annual monitoring is appropriate given the consistency of registry participation.

## Methods
Three anomaly detection approaches are implemented:


1. **Z‑Score Detection**: Year‑over‑year differences and percent changes are approximately normally distributed. Facilities are flagged when the absolute z‑score exceeds 1.5.


2. **Percent Change Thresholds**: Facilities are flagged when the absolute percent change between years exceeds 50%, reflecting operationally meaningful deviations.


3. **Negative Binomial–Poisson Prediction Interval**: A distribution‑based prediction interval is computed using a hybrid NB–Poisson model. NB is used when variance exceeds the mean. Poisson is used when variance is approximately equal to the mean. A facility is flagged when an observed count falls outside the predicted interval.

Quasi‑Poisson is not used because it is a quasi‑likelihood approach without a defined probability mass function or quantile function.

## Purpose

This project supports epidemiologic surveillance and quality review activities by identifying unusual reporting patterns using a mixture of statistical theory and operational pragmatism.