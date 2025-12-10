# spTimer-project-reproduction
STAT 5430 Final Project — Bayesian Spatio-Temporal Modeling Using spTimer
1. Overview
This project reproduces and extends the core simulation study from Bakar & Sahu (2015), focusing on Bayesian spatio-temporal modeling using the spTimer R package.
The original paper compares three models:
GP – Gaussian Process spatio-temporal model
AR – First-order autoregressive temporal model
GPP – Gaussian Predictive Process model (reduced-rank approximation)
This repository includes:
Full reproduction of the simulation study (12×12 and 55×55 spatial grids)
Full model refitting using spTimer (GP, AR, GPP)
Detailed performance comparisons (MAE, RMSE, R², coverage)
Three novel extension experiments that go beyond the original paper
This repository satisfies all requirements for the STAT 5430 Final Project:
Written report
Code repository with reproducibility
Final presentation slides
Substantial computational components
2. Extensions Introduced in This Project
In addition to reproducing the original simulation study, three computational extensions were conducted:
Extension 1 — Sampling Design (Grid vs Random)
We study how irregular spatial sampling affects prediction accuracy compared with regular grid sampling.
Random sampling increases variance in inter-site distances
Grid sampling achieves ≈5–8% lower MAE
Interpretation: Better spatial coverage improves kriging-style prediction
Extension 2 — Signal-to-Noise Ratio (SNR) Sensitivity
We vary SNR by adjusting process variance (sig2eta) and measurement error (sig2eps).
AR model is most sensitive to SNR; GPP remains stable
Extension 3 — Covariance Misspecification (Separable vs Nonseparable)
We simulate data from a nonseparable Gneiting kernel but fit a separable GP model.
Results show ≈5.5% loss in MAE
Coverage rates decline under misspecification
Demonstrates the importance of correctly capturing space–time interaction
4. How to Run the Analysis (you should change the address off source to your address)
1 data_simulation2 final.R
2 model_fitting2 final.R
3.nonseparable_cov3 final.R
4.covariance test2 final.R
5 grid_vs_random final.R
6 model_comparison2 final.R
7.Nonseparable Test with Robust Prediction Extraction2 final.R
8 original code from the paper and based on that I revised part of the code.R
5. Reproduction Results (Summary)
AR model performs worst across MAE/RMSE
GP achieves the best accuracy on small grids
GPP scales to 55×55 while retaining reasonable accuracy
6. Extension Experiment Results (Summary) 
The study compared two spatial sampling designs and two types of space–time correlation models.
A regular lattice provides full spatial coverage and near‑optimal spacing for Gaussian process interpolation on smooth fields.
A random sampling plan mimics realistic monitoring situations with clustered, irregular spacing.
The nonseparable correlation model shows strong diagonal dominance and abrupt transitions, capturing complex space–time interactions (e.g., pollutant diffusion, disease spread). The separable model yielded slightly lower correlation coefficients but more stable predictions.
Performance comparison showed:
MAE increased by about 5.5% (statistically significant).
RSE increased by about 1.4%.
R² remained nearly unchanged.
Overall, the results indicate that even similar overall correlation behavior can lead to different predictive accuracy, and the precise model specification substantially affects prediction reliability.
7. Documents
final paper
8. Experimental Scripts (Development History)
These files include the ones that I failed to achieve the successful outcomes







3.
