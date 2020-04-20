# datamoderate study

## Folder overview
Simulation study located in `simulation` folder with joint R code `simulations.R`.

### Folder `model_runs`:
* Life history scenarios (e.g. `shs`)
* Recruitment variability scenario (e.g. `LowSigmaR`)
* Files used to generate data with `ss3sim` (e.g. `cases`, `em`, `om`)
* Simulation scenario - will vary by life history, recruitment variability, and any alternate F scenario (e.g. `L1000-F1-I0-A0-E0-shs`)
* Simulation replicate, or iteration (e.g. `1`)
* Output files from `ss3sim` generated population -- true values found in `om`, ignore `em`
 * True population values in `om/Report.sso`
 * Length composition in `om/data.ss_new`
* Sample individuals per year to generate length composition (e.g. `N1000` for 1000 samples per year)
 * `perfect1000.ss` used `SS_splitdat` to find perfect information on the length composition
 * `ss3.dat` sampled N number of individuals from the perfect length composition using multinomial
 * Number of years of length data for estimation models -- all use `ss3.dat` but will adjust the number of years of length to include (e.g. `L100` uses 100 years, `L75` uses 75 years, `L1` uses last year only)
* Recruitment estimation scenarios, with all files for model run within. 

### Folder `R`
Includes some helper functions.

## Some additional details on the methods
Operating model uses R package `ss3sim`.
* Four life history scenarios e.g. `simulations/lifehistories`
 * shs = shorter-lived (max age = 30 years), slower-growing (expected to live at asymptotic length for final 10% of life)
 * shf = shorter-lived, faster-growing (expected to live at asymptotic length for final 50% of life)
 * los = longer-lived (max age = 60 years), slower-growing
 * lof = longer-lived, faster-growing
 * Each with variable M and k but sharing linf = 55 cm, length at 50% selectivity = 36.3 cm, h = 0.7, t0 = -1
* One fishing mortality time series, representative of U.S. West coast nearshore stocks. 
* Two recruitment scenarios: `LowSigmaR` sigmaR = 0.4 and `HighSigmaR` sigmaR = 0.8.
* 100 simulation replicates of each life history and recruitment scenario.
* Used to create the "true population", with values in `Report.sso` and information on length structure in `data.ss_new`

Data generation -- multinomial to sample from length structure 
* For each life history scenario and simulation replicate, generate length data from each year (e.g. `N1000`, `N1000`, `N50`):
 * 1000 samples (for checking)
 * 100 samples (more representative)
 * Changing from 100 to 50 samples over the data series
* Then choose the number of years to include in the model:
 * All 100 years of length data
 * Final 75 years, 20 years, 10 years, and 1 year, all subset from the same sampling procedure so that the length data is the same, only the number of years varies.

Estimation model
* Stock Synthesis (and will be tested with LIME)
* Recruitment estimation scenarios:
 * Unadjusted - using the default values for bias ramp from `ss3sim`
 * Bias adjusted - using `SS_fitbiasramp` to use estimated bias ramp parameters
 * No estimation - do no estimate recruitment deviates.

Performance
* Bias (median relative error) and precision (median absolute relative error)
* Interval coverage (proportion of iterations where true value lies within 50% confidence intervals;  nominal coverage would be equal to 50%. Scenarios greater than 50% tend to over-estimate uncertainty, while scenarios less than 50% tend to under-estimate uncertainty.)
 
