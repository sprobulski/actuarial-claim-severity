# Actuarial Modeling of Claim Severity and Heavy Tails

## Project Overview
This project focuses on statistical and actuarial modeling of individual claim severity using the `frecomfire` dataset. The analysis assumes that all claims are independent and identically distributed (i.i.d.), allowing for the unconditional modeling of severity through a single global distribution (without covariates). The primary goal is to accurately model heavy-tailed data and estimate extreme risk metrics for commercial fire losses. The workflow utilizes custom maximum likelihood estimations, the construction of complex spliced distributions, and rigorous goodness-of-fit testing tailored for extreme value theory.

## Documentation & Presentation
The repository includes a comprehensive presentation (`presentation/Claim_Severity_Presentation.pdf`) that provides:
* **Theoretical Background:** Mathematical formulations of heavy-tailed distributions, Extreme Value Theory (EVT), and risk measures (VaR, ES).
* **Implementation Details:** A breakdown of the Expectation-Maximization (EM) algorithm, parameter estimation, and threshold selection for spliced models.
* **Results Analysis:** A visual and quantitative comparison of classical distributions versus spliced models in capturing extreme tail events.

## Key Features & Methodology
* **Distributions Implemented:** Classical distributions (Gamma, Pareto, Fréchet, Log-logistic), Spliced models (Fréchet-Loglogistic, Fréchet-Pareto), and an Exponential distribution estimated via the EM algorithm.
* **Custom Implementations:** 
  * Manual log-likelihood functions and optimization routines (MLE and Moment Matching).
  * Custom architecture for spliced distributions, seamlessly connecting the body (Fréchet) and the heavy tail (Pareto/Log-logistic) at a predefined threshold (u = 40,000).
  * Analytical derivation and implementation of the Expectation-Maximization (EM) algorithm for missing/truncated tail data.
* **Actuarial Evaluation:** Tail behavior analysis (Hill Estimator, GPD), Value-at-Risk (VaR), Expected Shortfall (ES), Parametric and Non-parametric Bootstrap for confidence intervals, Probability Integral Transform (PIT), Exposure Curves, and Anderson-Darling/Kolmogorov-Smirnov statistical tests.

## Data
The dataset contains commercial fire insurance claims in France, sourced from the `CASdatasets` package (`frecomfire`). 
* **Target Variable:** `ClaimCost2007`, characterized by extreme asymmetry and a very heavy right tail.
* **Feature Engineering:** Logarithmic transformations and empirical tail smoothing via Gaussian Kernel Density Estimation (KDE).

## Key Results
The Spliced Fréchet-Pareto model significantly outperformed classical distributions, accurately balancing the body of the distribution and the extreme heavy tail (lowest AIC, best A-D statistic). Classical models like Gamma failed to capture extreme catastrophic events, while a pure Fréchet distribution overestimated the expected value. Furthermore, the bootstrap analysis revealed the extreme sensitivity of the Expected Shortfall (ES) metric to individual catastrophic claims within the historical sample.
