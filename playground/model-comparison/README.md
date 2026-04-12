# Model Comparison Playground

This folder contains the final expanded baseline working-model comparison used in
the thesis. All six methods use the upper-`tau` interval form `[lower(x), tau]`.

## Final Scripts

- `00_survreg_working_models.R`: Weibull and log-normal plug-in and conformal
  helpers.
- `01_run_model_comparison.R`: final six-method simulation driver.
- `02_analysis_model_comparison.qmd`: analysis/report script.
- `03_render_model_comparison.R`: HTML render helper.

Shared dependencies from `section-1-simulation/`:

- `00_utils.R`
- `01_dgp_weibull_ph.R`
- `02_conformal_cox_tau.R`
- `09_cox_plugin_Tstar.R`

## Final Results Kept

- `results/baseline-full-v5-all-upper/`
  - raw replicate outputs,
  - derived summary tables,
  - final PDF figures,
  - one convenience HTML report snapshot.

The tradeoff and direct-interval PDFs from this folder are the ones copied into
`ScM-thesis-writing/Figures/`.

## Reproduction

```bash
Rscript playground/model-comparison/01_run_model_comparison.R \
  --tag=baseline-full-v5-all-upper \
  --n_rep=30 \
  --n=2000 \
  --B=1500 \
  --seed0=20260321 \
  --tau_quantiles=0.70,0.80,0.90,0.99 \
  --wrong_dist=lognormal

Rscript playground/model-comparison/03_render_model_comparison.R \
  --results_dir=playground/model-comparison/results/baseline-full-v5-all-upper \
  --output_file=model-comparison-analysis.html
```
