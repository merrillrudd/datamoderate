# datamoderate study

## Folder overview
Simulation study located in `sim_new` folder with R code `sim_SS.R`.

### Model runs in life history scenario folders
* Life history scenarios (e.g. `short_slow`)
* Fishing mortality scenario (e.g `F1`)
* Recruitment variability scenario (e.g. `LowSigmaR`)
* Files used to generate data in folder `files`
* Simulation replicate, or iteration (e.g. `1`)
* Output files for generated population -- true values found in `om`
 * True population values in `om/Report.sso`
 * Length composition in `om/data.ss_new`
* Sample individuals per year to generate length composition
 * `perfect.ss` used `SS_splitdat` to find perfect information on the length composition
 * `ss3.dat` sampled N number of individuals from the perfect length composition using multinomial
 * Number of years of length data for estimation models -- all use `ss3.dat` but will adjust the number of years of length to include (e.g. `L100` uses 100 years, `L75` uses 75 years, `L1` uses last year only)
* Recruitment estimation scenarios, with all files for model run within. 

### Folder `R`
Includes some helper functions.

## Some additional details on the methods
Operating model runs initial life history, F, and recruitment scenario files without hessian to generate true population.
* Four life history scenarios
 * short_slow = shorter-lived (max age = 30 years), slower-to-Linf (expected to live at asymptotic length for final 10% of life)
 * short_fast = shorter-lived, faster-to-Linf (expected to live at asymptotic length for final 50% of life)
 * long_slow = longer-lived (max age = 60 years), slower-to-Linf
 * long_fast = longer-lived, faster-to-Linf
 * Each with variable M and k but sharing linf = 55 cm, length at 50% selectivity = 36.3 cm, h = 0.7, t0 = -1
* One fishing mortality time series, representative of U.S. West coast nearshore stocks. 
* Two recruitment scenarios: `LowSigmaR` sigmaR = 0.4 and `HighSigmaR` sigmaR = 0.8 (also explored deterministic)
* 100 simulation replicates of each life history, F, and recruitment scenario.
* Used to create the "true population", with values in `Report.sso` and information on length structure in `data.ss_new`

Data generation -- multinomial to sample from length structure 
* For each life history scenario and simulation replicate, generate length data from each year (e.g. `perfect`, `N000`, `N50`):
 * perfect information with 1000 length samples
 * 100 samples (more representative)
 * Changing from 100 to 50 samples over the data series
* Then choose the number of years to include in the model:
 * All 100 years of length data with perfect information
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
 
