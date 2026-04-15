# Challenging DGP Comparison

This playground compares the same six upper-`tau` methods used in the thesis
model-comparison section, but under harder finite-window DGPs where the Cox
working model is misspecified.

Current scenarios:

- `aft_lognormal_linear`: linear log-normal AFT DGP; `Log-Normal` is the correct
  parametric working model, while `Cox` and `Weibull` are misspecified.
- `aft_lognormal_nonlinear`: nonlinear log-normal AFT DGP; all three working
  models are misspecified.

Core files:

- `00_dgp_challenging.R`: challenging DGP definitions.
- `01_run_challenging_model_comparison.R`: simulation driver.
- `02_analysis_challenging_model_comparison.qmd`: summary plots and tables.
- `03_render_challenging_model_comparison.R`: HTML renderer.

Dependencies loaded from the main project:

- `section-1-simulation/00_utils.R`
- `section-1-simulation/01_dgp_weibull_ph.R`
- `section-1-simulation/02_conformal_cox_tau.R`
- `section-1-simulation/09_cox_plugin_Tstar.R`
- `playground/model-comparison/00_survreg_working_models.R`
