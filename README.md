# 'Polarization on Social Media' replication package

Replication materials for **Polarization on Social Media: Micro-Level Evidence and Macro-Level Implications**, JASSS 27(1), 7 (2024). [Published paper](https://www.jasss.org/27/1/7.html), [DOI](https://doi.org/10.18564/jasss.5298).

This compact package starts from the project's cleaned, model-ready observations. It contains de-identified analysis data, the original statistical model functions, a portable analysis pipeline, selected historical ABM code, and a standalone ABM implementation. Original project files have not been changed. Participant anonymity is the aim; the public paper remains cited for attribution.

## Quick start

Open `master.R` in RStudio. Set the working directory to this folder, install the three plotting dependencies once, then work through the clearly marked sections line by line:

```r
install.packages(c("ggplot2", "patchwork", "ggExtra"))
```

For a complete batch run from this folder:

```sh
Rscript master.R
```

The master script displays plots in interactive sessions and saves PNG/PDF versions in `output/plots/`. It contains the plotting code directly, keeps named plot objects available for editing, and explains the data filters, parameters and figure mappings. It adds distributions of opinions, ideology and signed opinion change, pre/post comparisons, and follow-up descriptives. Missing simulation CSVs are generated using Python 3; set `run_missing_simulations <- FALSE` to run only the available plots. Existing simulation CSVs are reused.

The original posterior draws are in `data/posterior/`. They are analysis inputs required for the predictions and calibrated ABM, not validation artifacts. No validation reports or comparison files are distributed. The default workflow summarizes these draws rather than refitting. All generated files go in `output/` and can be deleted and recreated.

The separate table/statistics entry point remains available:

```sh
Rscript code/analyze.R
```

This calculates descriptive/posterior summaries and WAIC. For LOO, install the R package `loo` first. The original `code/influence_functions.R` and `abm/plot.R` scripts remain available, but the master script does not depend on them.

## Install defSim for the notebook

The first code cell of `abm/original_experiments.ipynb` installs defSim from [its GitHub repository](https://github.com/defSim/defSim) into the active notebook kernel:

```python
%pip install "git+https://github.com/defSim/defSim.git" matplotlib seaborn
```

Alternatively, install it from a terminal into the Python environment used by Jupyter:

```sh
python3 -m pip install "git+https://github.com/defSim/defSim.git" matplotlib seaborn notebook
python3 -m notebook abm/original_experiments.ipynb
```

Git and internet access are needed. Restart the notebook kernel after installation if prompted. This follows the repository's GitHub installation method and installs its current version; it does not reconstruct the unrecorded historical defSim environment. The standalone `abm/run.py` used by the master R script does not require defSim.

## Contents

| Location | Purpose |
| --- | --- |
| `data/study1.csv` | 312 topic-level observations from the cleaned Study 1 data, including follow-up opinions where observed. |
| `data/study2.csv` | 474 topic-level observations from the cleaned Study 2 data, including controls. |
| `data/README.md` | Variable definitions, missing values, exclusions, and anonymization details. |
| `master.R` | Annotated, interactive plotting workflow: paper figures, additional descriptives, and optional simulation generation. |
| `code/common.R` | Data loading, original model-specific sample filters, predictors, and reconstruction of pointwise log likelihoods. |
| `code/models_study1.R`, `code/models_study2.R` | Original JAGS model functions, extracted from the project's `04_models` scripts. Statistical specifications, priors, chain seeds, and iteration defaults are preserved. |
| `code/analyze.R` | Posterior summaries, WAIC, optional LOO and diagnostics, descriptive statistics, and simplified versions of empirical Figures 3–5. |
| `code/refit.R` | Refit the 12 models underlying the main and follow-up tables using the de-identified data. |
| `code/influence_functions.R` | Original influence-function plotting code, with output paths made relative. |
| `abm/run.py` | Portable one-topic ABM with explicit seeds and alternative documented experiment designs. |
| `abm/plot.R` | Plot simulation means using base R. |
| `abm/original_experiments.ipynb` | Selected original model definitions and experiments for provenance; outputs, execution counts, metadata, and personal paths removed. Requires the historical defSim environment; it is **not** the recommended run-all entry point. |
| `data/posterior/*.csv.gz` | Original posterior parameter draws: 3 chains × 2,000 retained draws per model; parameters, chain numbers and retained-draw indices only. |

## Refit from the released observations

Install **JAGS 4** separately and install the R package `runjags`:

```r
install.packages("runjags")
```

Then, from the package root:

```sh
Rscript code/refit.R
Rscript code/analyze.R output/refit
python3 abm/run.py --profile calibrated --draws output/refit
```

Each original model function requests 3 chains, 4,000 burn-in iterations, 2,000 retained samples per chain and thinning of 5. Initial RNG types and seeds are specified inside the functions. Additional adaptation is controlled by runjags/JAGS. Refitting may take substantially longer than summarizing the archived draws. Differences across software versions, adaptation, and the shuffled release-data order mean refits need not reproduce the original chains exactly.

## Mapping to the paper

| Paper result | Replication route |
| --- | --- |
| Table 1 | `output/descriptive_statistics.csv`; uses the original descriptive complete-case filter. |
| Table 2 | `m10`, `m11`, `m12` in `posterior_summary.csv` and `model_fit.csv`; N = 311. |
| Table 3 | `m10.lag`, `m11.lag`, `m12.lag`; N = 124. |
| Table 4 | `m20`, `m21`, `m22`; N = 306. |
| Table 5 (appendix reference models) | `m10.ref`, `m11.ref`, `m12.ref`; N = 124. |
| Figures 3 and 5 | `master.R` section 4: absolute opinion change versus perceived distance, with loess curves and marginal histograms. |
| Figure 4 | `master.R` section 5: response functions from Model 2 posterior means, at ideological endpoints and both moral frames. |
| Figure 1 | ABM `illustration` profile: 16 homophily values × 2 bubble sizes × 10 repetitions = 320 runs, 1,000 interactions each. New seeded realizations, not the original image. |
| Figure 6 | See the design comments in `master.R` section 7. `calibrated` preserves the notebook's stochastic-calibration experiment; `caption` provides an explicitly reconstructed caption-based design. Neither is claimed to recover the exact historical Figure 6 runs. |
| Theoretical influence-function illustrations | `master.R` section 6 generates Figures 2 and 8 from the local plotting formulas, plus a separately labeled version of published equation 13, which differs from the original nonlinear plotting code. These are not pixel-identical reconstructions. |

The original R coefficient named `gamma` corresponds to the paper's estimated distance coefficient **beta**, not the ABM homophily ratio. ABM homophily is beta/alpha. `zeta` is the Study 2 ideological-distance interaction. In Study 1, `b_moral`, `b_right`, and `b_moralright` are the three added interactions.

After refitting, set `draws_dir <- "output/refit"` near the top of `master.R`. Delete existing ABM CSVs before rerunning if the simulations should use the new calibration. The script labels the notebook and caption-based Figure 6 designs separately and reports which coverage calculation is plotted. Figure 7 is a study-design schematic and is not regenerated. The package begins with cleaned data, so raw-survey processing and participant-linked analyses are outside its scope.

## AI statement

This replication package was created with the assistance of OpenAI's ChatGPT, which helped organizing the original code (no AI used in creating that code) and draft and edit the README and other documentation. The authors take full responsibility for the content and any errors in the package.
