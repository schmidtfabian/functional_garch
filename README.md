# On Functional GARCH Models and Their Application to Intraday Value-at-Risk
This repository contains the codes of my master thesis "On Functional GARCH Models and Their Application to Intraday Value-at-Risk".

Quick overview over the files contained in this repository:
- main_simulation_1.R and main_simulation_2.R contain the simulation codes that produce the results of the simulation section.
- main.tex contains the Latex code that produced the final pdf version of the thesis.
- bibliography.bib contains the references.
- custom_chicago.bst contains the template of a custom bibliography style based on the chicago.bst template.
- Master_Thesis_Fabian_Schmidt.pdf is the pdf version of the thesis.
- The subfolder /R contains all R functions used in the simulations.

Here is an overview over the files contained in the /R folder:
- baseline_functions.R, FPC.R and PCA_compute compute the baseline functions, basis functions and Principal components used for the FGARCH QML, FGARCH LS and PC-GARCH model respectively.
- simulate_fgarch_1.R and simulate_fgarch_2.R contains the R functions that simulate the burn-in, training and evaluation sample in the Monte Carlo Simulations.
- fgarch_est_ls.R, fgarch_est_qml.R and pc_garch_estimate contain the R functions for the estimation of the parameters of the FGARCH LS, FGARCH QML and PC-GARCH model respectively.
- fgarch_est_ls_gradient.R is an alternative LS estimation using gradients and Jacobians that was NOT used to obtain the results found in the simulation section of the thesis.
- epsilon_functions.R contain three functions meant to simulate the innovation functions of the FGARCH(1,1) process. Two of those functions were not used to obtain the results found in the simulation section of the thesis.
- calculate_fitted_values.R calculates the fitted values for all three models.
- compute_test_statistic.R computes the mean p-values and standard deviation of the p-values and conditionally the root squared mean deviations from the true functional parameters.
- calculate_bootstrap.R calculates the bootstrap functions from the fitted innovation functions provided.
- forecast_VAR.R and forecast_pc_garch forecasts the value-at-risk values for the FGARCH model and the PC-GARCH model respectively.