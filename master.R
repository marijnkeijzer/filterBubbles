# Polarization on Social Media: guided replication -------------------------------
# Paper: https://www.jasss.org/27/1/7.html
# Open this file in RStudio and run sections / individual lines in order.
# Or: Rscript master.R, from the replication-package directory.
# All plotting code is here. Named plot objects stay available for inspection.
# Plots display interactively and are saved as PNG and PDF in output/plots.
# Figure numbers identify analyses, not pixel-identical copies of the artwork.
# Figure 7 is a study-design schematic, not an analysis output.

# 1. Setup: edit these settings first -------------------------------------------
# RStudio: Session > Set Working Directory > To Source File Location.
# Alternatively replace "." with the path to this replication-package folder:
package_dir <- "."
setwd(package_dir)
stopifnot(file.exists("data/study1.csv"), file.exists("data/study2.csv"))

# One-time installation (uncomment if necessary):
# install.packages(c("ggplot2", "patchwork", "ggExtra"))
required <- c("ggplot2", "patchwork", "ggExtra")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "))
library(ggplot2)
library(patchwork)

# Original posterior draws; change to "output/refit" after fitting new models.
draws_dir <- "data/posterior"
# TRUE generates missing ABM CSVs using the included standard-library Python model.
# FALSE skips missing simulations; empirical and theoretical plots still run.
# Existing CSVs are reused. Delete them to change seed or calibration draws.
run_missing_simulations <- TRUE
simulation_seed <- 20261117L
python <- Sys.which("python3")
if (!nzchar(python)) python <- Sys.which("python")
plot_dir <- "output/plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.minor = element_blank(), legend.position = "bottom"))

# In an interactive session, typing any plot object's name displays it again.
# Batch runs avoid opening an extra default device / creating Rplots.pdf.
save_plot <- function(plot, name, width = 8, height = 5) {
  if (interactive()) print(plot)
  ggsave(file.path(plot_dir, paste0(name, ".png")), plot, width = width,
         height = height, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, paste0(name, ".pdf")), plot, width = width,
         height = height, bg = "white")
}

# 2. Read and inspect the data --------------------------------------------------
# Rows are participant-topic observations, not participant counts.
# Identifiers were removed; matching respondents across topics is not possible.
# t0 = before stimulus; t1 = immediately after; t2 = Study 1 follow-up.
# Opinions and sti_lr are on 0-1 scales; leftright is on a 0-100 scale.
study1 <- read.csv("data/study1.csv", na.strings = "NA")
study2 <- read.csv("data/study2.csv", na.strings = "NA")
str(study1)
str(study2)
table(study1$topic)
table(study2$topic)

observations <- rbind(
  data.frame(study = "Study 1", study1[, names(study2)]),
  data.frame(study = "Study 2", study2)
)
observations$topic <- factor(observations$topic, levels = c("dev", "tax"),
                             labels = c("Development aid", "Taxation"))
observations$distance <- abs(observations$sti_auth - observations$t0_gen)
observations$change <- observations$t1_gen - observations$t0_gen
observations$abs_change <- abs(observations$change)
observations$boundary <- factor(ifelse(observations$t1_gen %in% c(0, 1),
                                       "At boundary", "Within scale"))

# 3. Additional descriptives: distributions and observed change ------------------
# These complement the paper rather than recreating numbered figures.
# Available observations are used; missing values are dropped plot by plot.
opinion_long <- rbind(
  data.frame(observations[c("study", "topic")], time = "Before", opinion = observations$t0_gen),
  data.frame(observations[c("study", "topic")], time = "After", opinion = observations$t1_gen)
)
opinion_long$time <- factor(opinion_long$time, levels = c("Before", "After"))
p_opinions <- ggplot(opinion_long, aes(opinion, fill = time)) +
  geom_histogram(binwidth = .1, boundary = 0, position = "identity", alpha = .45,
                 na.rm = TRUE) + facet_grid(study ~ topic) +
  scale_fill_grey(start = .2, end = .8) +
  labs(x = "Opinion (0-1)", y = "Topic-level observations", fill = "Measurement")
save_plot(p_opinions, "descriptive_opinion_distributions")

p_pre_post <- ggplot(observations, aes(t0_gen, t1_gen)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
  geom_point(alpha = .4, size = 1, na.rm = TRUE) + facet_grid(study ~ topic) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Opinion before stimulus", y = "Opinion after stimulus")
save_plot(p_pre_post, "descriptive_pre_post", height = 7)

p_change <- ggplot(observations, aes(change)) +
  geom_histogram(binwidth = .05, boundary = 0, fill = "grey75", colour = "white", na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = 2) + facet_grid(study ~ topic) +
  labs(x = "Signed opinion change (after - before)", y = "Topic-level observations")
save_plot(p_change, "descriptive_signed_change")

p_ideology <- ggplot(observations, aes(leftright)) +
  geom_histogram(binwidth = 10, boundary = 0, fill = "grey65", colour = "white", na.rm = TRUE) +
  facet_grid(study ~ topic) + labs(x = "Left-right identification (0-100)",
                                 y = "Topic-level observations")
save_plot(p_ideology, "descriptive_ideology")

followup <- study1[complete.cases(study1[c("t0_gen", "t1_gen", "t2_gen")]), ]
followup_long <- rbind(
  data.frame(topic = followup$topic, wave = "Immediate", change = followup$t1_gen - followup$t0_gen),
  data.frame(topic = followup$topic, wave = "Follow-up", change = followup$t2_gen - followup$t0_gen)
)
p_followup <- ggplot(followup_long, aes(change, fill = wave)) +
  geom_histogram(binwidth = .05, boundary = 0, position = "identity", alpha = .45) +
  facet_wrap(~topic) + scale_fill_grey(start = .2, end = .8) +
  labs(x = "Opinion change from baseline (complete follow-up subset)",
       y = "Topic-level observations", fill = "Measurement")
save_plot(p_followup, "descriptive_followup")

# 4. Figures 3 and 5: opinion change versus perceived distance --------------------
# Original figures use available model-ready observations, including Study 2
# controls. This is not the restricted model-fitting sample.
# Loess is descriptive (span .75), not the Bayesian influence function.
# Marginal histograms show the distribution along each axis, as in the paper.
plot_observed_change <- function(study_label) {
  d <- observations[observations$study == study_label &
                      complete.cases(observations[c("distance", "abs_change")]), ]
  scatter <- ggplot(d, aes(distance, abs_change)) +
    geom_smooth(method = "loess", formula = y ~ x, span = .75, se = TRUE,
                colour = "black", fill = "grey80", linewidth = .6) +
    geom_point(aes(shape = boundary), size = 1.5) +
    scale_shape_manual(values = c("At boundary" = 2, "Within scale" = 1)) +
    labs(x = "Perceived distance to stimulus", y = "Absolute opinion change",
         shape = "Post-stimulus opinion")
  ggExtra::ggMarginal(scatter, type = "histogram", fill = "grey85", colour = "white")
}
figure3 <- plot_observed_change("Study 1")
save_plot(figure3, "figure3_study1_observed_change", width = 7, height = 6)
figure5 <- plot_observed_change("Study 2")
save_plot(figure5, "figure5_study2_observed_change", width = 7, height = 6)

# 5. Figure 4: response functions by ideology and moral alignment -----------------
# Posterior-mean model predictions, not empirical subgroup averages.
# R coefficient `gamma` is the paper's estimated beta (distance coefficient).
# Author opinion = 0; ideology evaluated at endpoints 0 (left) and 1 (right).
read_parameters <- function(model) {
  file <- file.path(draws_dir, paste0(model, ".csv.gz"))
  if (!file.exists(file)) stop("Missing parameter draws: ", file)
  draws <- read.csv(gzfile(file))
  draws[setdiff(names(draws), c("chain", "iteration"))]
}
m12_draws <- read_parameters("m12")
m12_mean <- vapply(m12_draws, mean, numeric(1))
m12_mean

prediction_grid <- expand.grid(before = seq(0, 1, length.out = 501), moral = 0:1, ideology = 0:1)
prediction_grid$effective_alpha <- with(prediction_grid,
  m12_mean["alpha"] + m12_mean["b_moral"] * moral +
    m12_mean["b_right"] * ideology + m12_mean["b_moralright"] * moral * ideology)
prediction_grid$delta <- -prediction_grid$before
prediction_grid$after <- with(prediction_grid,
  before + effective_alpha * delta - m12_mean["gamma"] * abs(delta) * delta)
prediction_grid$after <- pmin(1, pmax(0, prediction_grid$after))
prediction_grid$alignment <- ifelse(prediction_grid$moral == prediction_grid$ideology, "Aligned", "Disaligned")
prediction_grid$orientation <- ifelse(prediction_grid$ideology == 0, "Left", "Right")
figure4 <- ggplot(prediction_grid, aes(before, after, colour = alignment, linetype = orientation)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey75") + geom_line(linewidth = .8) +
  scale_colour_grey(start = .1, end = .6) + coord_equal() +
  labs(x = "Pre-stimulus opinion", y = "Predicted post-stimulus opinion",
       colour = "Moral foundation", linetype = "Ideological identification") +
  theme(legend.box = "vertical")
save_plot(figure4, "figure4_moral_alignment", width = 7, height = 6)

# 6. Figures 2 and 8: theoretical weight and response functions --------------------
# Alpha = 1 and author opinion = 0, as in the figure captions.
# Figure 2: weight = 1 - gamma * distance; gamma = .5, 1, 2, 4.
# Figure 8: local plotting code uses (1-distance)^beta, but published equation 13
# prints 1-distance^beta. They differ. Preserve the local plotting convention
# for figure8 and show the printed-equation alternative in a separate plot.
theory <- expand.grid(distance = seq(0, 1, length.out = 501), parameter = c(.5, 1, 2, 4))
theory$curve <- factor(theory$parameter)
theory$linear_weight <- with(theory, 1 - parameter * distance)
theory$linear_after <- with(theory, pmin(1, pmax(0, distance - linear_weight * distance)))
theory$nonlinear_weight <- with(theory, (1 - distance)^parameter)
theory$nonlinear_after <- with(theory, distance - nonlinear_weight * distance)
theory$equation13_weight <- with(theory, 1 - distance^parameter)
theory$equation13_after <- with(theory, distance - equation13_weight * distance)

# .data[[...]] selects a named data column; it does not change the data.
plot_theory <- function(weight_column, response_column, legend_title) {
  weight_plot <- ggplot(theory, aes(distance, .data[[weight_column]], linetype = curve)) +
    geom_hline(yintercept = 0, colour = "grey75") + geom_line(linewidth = .7) +
    labs(x = "Opinion distance", y = "Influence weight", linetype = legend_title)
  response_plot <- ggplot(theory, aes(distance, .data[[response_column]], linetype = curve)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey75") + geom_line(linewidth = .7) +
    labs(x = "Pre-stimulus opinion", y = "Post-stimulus opinion", linetype = legend_title)
  (weight_plot | response_plot) + plot_layout(guides = "collect")
}
figure2 <- plot_theory("linear_weight", "linear_after", "gamma")
save_plot(figure2, "figure2_linear_influence", width = 10, height = 4)
figure8 <- plot_theory("nonlinear_weight", "nonlinear_after", "beta")
save_plot(figure8, "figure8_nonlinear_original_code", width = 10, height = 4)
p_equation13 <- plot_theory("equation13_weight", "equation13_after", "beta")
save_plot(p_equation13, "supplement_equation13_convention", width = 10, height = 4)

# 7. Figures 1 and 6: run/read simulations and plot their outcomes -----------------
# The portable Python model does not require defSim. See the notebook's first
# cell or README to install defSim for the historical notebook instead.
# illustration: 320 runs, N=50, 1,000 interactions, original parameter grid.
# calibrated: notebook design, 40 runs, N=50, 5,000 interactions.
# caption: reconstructed design, 10 repeats x 1,000 interactions per bubble size
# 0-50. Paper text instead mentions 10,000 runs. Figure 6 alternatives are labeled;
# neither claims to recover the exact historical simulation sample.
profiles <- c("illustration", "calibrated", "caption")
for (profile in profiles) {
  csv_file <- file.path("output", paste0("abm_", profile, ".csv"))
  if (!file.exists(csv_file) && run_missing_simulations) {
    if (!nzchar(python)) stop("Python 3 not found. Install it or set run_missing_simulations <- FALSE.")
    status <- system2(python, c(shQuote("abm/run.py"), "--profile", profile,
                               "--seed", simulation_seed, "--draws", shQuote(draws_dir)))
    if (status != 0) stop("Simulation failed for profile: ", profile)
  }
}

# Curves = means; ribbons = +/- 1.96 standard errors across independent runs.
# These normal-approximation ribbons differ from seaborn's original bootstrap.
# Coverage legacy retains the original reporter's indexing error; calibrated
# panels explicitly use the corrected interval-union measure instead.
plot_simulations <- function(profile) {
  file <- file.path("output", paste0("abm_", profile, ".csv"))
  if (!file.exists(file)) {
    message("Skipping ", profile, ": generate ", file, " first.")
    return(NULL)
  }
  sim <- read.csv(file)
  illustration <- profile == "illustration"
  measures <- c("spread", "dispersion", if (illustration) "coverage_legacy" else "coverage_interval")
  long <- do.call(rbind, lapply(measures, function(measure) {
    data.frame(x = if (illustration) sim$gamma else sim$bubble,
               group = if (illustration) paste("Bubble", sim$bubble) else "All runs",
               measure = measure, value = sim[[measure]])
  }))
  grouped <- split(long, interaction(long$x, long$group, long$measure, drop = TRUE))
  means <- do.call(rbind, lapply(grouped, function(z) {
    data.frame(x = z$x[1], group = z$group[1], measure = z$measure[1],
               mean = mean(z$value), se = sd(z$value) / sqrt(nrow(z)), n = nrow(z))
  }))
  means$measure <- factor(means$measure, levels = measures,
    labels = c("Spread", "Dispersion", if (illustration) "Coverage (legacy)" else "Coverage (interval union)"))
  p <- ggplot(means, aes(x, mean, linetype = group, group = group)) +
    geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                fill = "grey70", alpha = .3, linetype = 0) +
    geom_line(linewidth = .7) + facet_wrap(~measure, nrow = 1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = if (illustration) "Homophily (gamma)" else "Bubble size (agents)",
         y = "Mean polarization", linetype = NULL, subtitle = paste("Simulation design:", profile))
  if (!illustration) p <- p + guides(linetype = "none")
  p
}
figure1 <- plot_simulations("illustration")
if (!is.null(figure1)) save_plot(figure1, "figure1_simulation", width = 11, height = 4)
figure6_notebook <- plot_simulations("calibrated")
if (!is.null(figure6_notebook)) save_plot(figure6_notebook, "figure6_notebook_design", width = 11, height = 4)
figure6_caption <- plot_simulations("caption")
if (!is.null(figure6_caption)) save_plot(figure6_caption, "figure6_caption_design", width = 11, height = 4)

# 8. Optional: tables and fresh estimation ---------------------------------------
# Plotting is complete. Commands below are commented because fitting needs
# JAGS 4 / runjags and takes longer. Summaries do not need JAGS.
# For posterior/descriptive tables, WAIC and optional LOO:
# system2(file.path(R.home("bin"), "Rscript"), c("code/analyze.R", shQuote(draws_dir)))
# For LOO first install.packages("loo"). For fresh estimates, install JAGS 4 and:
# install.packages("runjags")
# source("code/refit.R")
# Restart this guide with draws_dir <- "output/refit" and delete existing ABM
# CSVs if simulations should use the new calibration. Original truncated
# likelihood and signed-sigma prior are preserved in code/models_study*.R.

message("Finished. Plots are in ", normalizePath(plot_dir), ".")
