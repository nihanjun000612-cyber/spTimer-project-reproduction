STAT 5430 Final Project — Bayesian Spatio-Temporal Modeling Using spTimer
1. Overview
   This project reproduces and extends the core simulation study from Bakar & Sahu (2015), focusing on Bayesian spatio‑temporal modeling using the spTimer R package.
Original Models Compared
GP — Gaussian Process spatio‑temporal model
AR — First‑order autoregressive temporal model
GPP — Gaussian Predictive Process model (reduced‑rank approximation)
This Repository Includes
Full reproduction of the simulation study (12×12 and 55×55 spatial grids)
Full model refitting using spTimer (GP, AR, GPP)
Detailed performance comparisons (MAE, RMSE, R², coverage)
Three computational extensions beyond the original paper
This repository satisfies the full requirements for the STAT 5430 Final Project:
Written report
Reproducible code repository
Final presentation slides
Substantial computational components

2. Extensions Introduced in This Project
Extension 1 — Sampling Design (Grid vs Random)
Examines how irregular spatial sampling affects prediction accuracy compared to regular grid sampling.
Random sampling increases variance in inter‑site distances.
Grid sampling achieves ≈5–8% lower MAE.
Interpretation: Better spatial coverage improves kriging‑style prediction.
Extension 2 — Signal‑to‑Noise Ratio (SNR) Sensitivity
Varies SNR by adjusting process variance (sig2eta) and measurement error (sig2eps).
AR model is most sensitive to SNR, while GPP remains more stable.
Extension 3 — Covariance Misspecification (Separable vs Nonseparable)
Simulates data from a nonseparable Gneiting kernel but fits a separable GP model.
Results show ≈5.5% loss in MAE and a small 1.4% increase in RSE.
3. How to Run the Analysis
(Change file paths to your local environment before running.)
1.data simulation2 final
2.model fitting2 final.R
3.nonseparable cov3 final.R
4.covariance test2 final.R
5.grid vs random final.R
6.Nonseparable Test with Robust Prediction Extraction2 final.R
7.Nonseparable_Test_with_Robust_Prediction_Extraction2_final.R
8.original code from paper revised.R
4. Reproduction Results (Summary)
Model	Performance Summary
AR	Performs worst across MAE/RMSE
GP	Best accuracy on smaller grids
GPP	Scales to 55×55 while retaining reasonable accuracy
5. Extension Experiment Results (Summary)
Two spatial sampling designs: regular lattice vs random sampling.
Two space–time correlation structures: separable and nonseparable models.
Findings:
Nonseparable model shows strong diagonal dominance and abrupt transitions — it better captures complex space–time interactions (e.g., pollutant diffusion, disease spread).
Separable model has slightly lower correlation coefficients but yields more stable predictions.
MAE increased by ~5.5% (statistically significant).
RSE increased by ~1.4%.
R² remained nearly unchanged.
Even when overall correlation behavior looks similar, correct model specification is crucial for improving predictive accuracy and model stability.
6. Documents
Bayesian Spatio-Temporal Modeling Using spTimer A Reproduction and Extension Study.docx/ original paper.pdf
7. Experimental Scripts (Development History)
Includes earlier or incomplete versions used during the development process.
