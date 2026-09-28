## -----------------------------------------------------------------------------
#| echo: false
#| cache: false
library(rpact)


## -----------------------------------------------------------------------------
#| message: false
library(crmPack)

dose_grid <- c(1, 5, seq(from = 10, to = 200, by = 1))


## -----------------------------------------------------------------------------
empty_data <- Data(doseGrid = dose_grid)


## -----------------------------------------------------------------------------
data <- empty_data |>
    update(x = 10, y = 0) |>
    update(x = 20, y = 1)


## -----------------------------------------------------------------------------
plot(data)


## -----------------------------------------------------------------------------
get_cov <- function(sigma0, sigma1, rho) {
    matrix(c(sigma0^2, rho * sigma0 * sigma1, rho * sigma0 * sigma1, sigma1^2), nrow = 2)
}
model <- LogisticLogNormal(
    mean = c(-3, 0.5),
    cov = get_cov(sigma0 = 1, sigma1 = 2, rho = 0.5),
    ref_dose = 20
)


## -----------------------------------------------------------------------------
model


## -----------------------------------------------------------------------------
mcmc_options <- McmcOptions(
    burnin = 1000,
    step = 10,
    samples = 1000,
    rng_kind = "Mersenne-Twister",
    rng_seed = 12345
)


## -----------------------------------------------------------------------------
prior_samples <- mcmc(model, data = empty_data, options = mcmc_options)
posterior_samples <- mcmc(model, data = data, options = mcmc_options)


## -----------------------------------------------------------------------------
names(posterior_samples@data)
alpha0_samples <- get(posterior_samples, "alpha0")
library(ggmcmc)
ggs_traceplot(alpha0_samples)
ggs_autocorrelation(alpha0_samples)


## -----------------------------------------------------------------------------
plot(prior_samples, model, empty_data)


## -----------------------------------------------------------------------------
plot(posterior_samples, model, data)


## -----------------------------------------------------------------------------
increments <- IncrementsRelativeDLT(
    intervals = c(0, 1),
    increments = c(1, 0.5)
)
increments


## -----------------------------------------------------------------------------
next_best <- NextBestNCRMLoss(
    target = c(0.2, 0.4), # loss: 0
    overdose = c(0.4, 0.7), # loss: 1
    unacceptable = c(0.7, 1), # loss: 2
    max_overdose_prob = 0.25,
    losses = c(1, 0, 1, 2) # underdosing has loss 1
)
next_best


## -----------------------------------------------------------------------------
cohort_size <- CohortSizeConst(3)


## -----------------------------------------------------------------------------
stop1 <- StoppingMinPatients(nPatients = 30)


## -----------------------------------------------------------------------------
stop2 <- StoppingTargetProb(target = c(0.2, 0.3), prob = 0.7)
stop3 <- StoppingMissingDose()
stopping <- stop1 | stop2 | stop3
stopping


## -----------------------------------------------------------------------------
max_dose <- maxDose(increments, data)
max_dose


## -----------------------------------------------------------------------------
next_best_dose <- nextBest(
    next_best,
    doselimit = max_dose,
    samples = posterior_samples,
    model = model,
    data = data
)


## -----------------------------------------------------------------------------
print(next_best_dose$plot_joint)


## -----------------------------------------------------------------------------
next_best_dose$value


## -----------------------------------------------------------------------------
fitted_df <- fit(posterior_samples, model, data)
loss_df <- next_best_dose$probs |>
    as.data.frame() |>
    dplyr::select(dose, posterior_loss)
result_df <- dplyr::left_join(fitted_df, loss_df, by = "dose")
knitr::kable(head(result_df), digits = 3)


## -----------------------------------------------------------------------------
#| results: false
design <- Design(
    data = empty_data,
    model = model,
    nextBest = next_best,
    increments = increments,
    cohort_size = cohort_size,
    stopping = stopping,
    startingDose = 5
)
design

## -----------------------------------------------------------------------------
#| echo: false
#| results: asis
design_text <- as.character(knitr::knit_print(design))
design_text <- gsub("(?m)^##\\s+", "### ", design_text, perl = TRUE)
design_text <- gsub("(?m)^###\\s+", "#### ", design_text, perl = TRUE)
cat(design_text)


## -----------------------------------------------------------------------------
examine(design)


## -----------------------------------------------------------------------------
model_adjusted <- LogisticLogNormal(
    mean = c(-3, -1),
    cov = get_cov(sigma0 = 4, sigma1 = 1, rho = 0),
    ref_dose = 100
)
prior_samples_adjusted <- mcmc(model_adjusted, data = empty_data, options = mcmc_options)
plot(prior_samples_adjusted, model_adjusted, empty_data)

design_adjusted <- design
design_adjusted@model <- model_adjusted
examine(design_adjusted)


## -----------------------------------------------------------------------------
scenario1 <- function(x) {
    plogis(-2 + 0.8 * log(x / 20))
}


## -----------------------------------------------------------------------------
curve(scenario1, from = 1, to = 200, ylim = c(0, 1))


## -----------------------------------------------------------------------------
sims <- simulate(
    design_adjusted,
    nsim = 100,
    truth = scenario1,
    mcmcOptions = mcmc_options,
    seed = 123,
    parallel = TRUE
)


## -----------------------------------------------------------------------------
backfill_simple <- Backfill(
    cohort_size = CohortSizeConst(3),
    max_size = 12,
    opening = OpeningMinCohorts(min_cohorts = 1),
    recruitment = RecruitmentUnlimited(),
    priority = "lowest"
)
backfill_simple


## -----------------------------------------------------------------------------
design_backfill <- design_adjusted
design_backfill@backfill <- backfill_simple


## -----------------------------------------------------------------------------
sims_summary <- summary(sims, truth = scenario1)
sims_summary


## -----------------------------------------------------------------------------
plot(sims)


## -----------------------------------------------------------------------------
plot(sims_summary)


## -----------------------------------------------------------------------------
#| message: false
coarseGrid <- c(25, 100, 300)
min_result <- MinimalInformative(
    dosegrid = coarseGrid,
    refDose = 100,
    logNormal = TRUE, # use a log-normal distribution for the slope parameter
    threshmin = 0.1, # the quantile for the low dose
    threshmax = 0.2, # the quantile for the high dose
    seed = 432, # for reproducibility
    control = list(max.time = 30) # limit the optimization time here
)


## -----------------------------------------------------------------------------
#| output-location: slide
matplot(
    x = coarseGrid,
    y = min_result$required,
    type = "o",
    lty = 1,
    col = 1
)
matlines(
    x = coarseGrid,
    y = min_result$quantiles,
    type = "o",
    lty = 2,
    col = 1
)
legend(
    "topleft",
    legend = c("Target quantiles", "Obtained quantiles"),
    lty = 1:2,
    col = 1
)


## -----------------------------------------------------------------------------
min_result$model


## -----------------------------------------------------------------------------
comp1 <- ModelParamsNormal(
    mean = c(-0.85, 1),
    cov = matrix(c(1, -0.5, -0.5, 1), nrow = 2)
)
comp2 <- ModelParamsNormal(
    mean = c(1, 1.5),
    cov = matrix(c(1.2, -0.45, -0.45, 0.6), nrow = 2)
)


## -----------------------------------------------------------------------------
normal_mix_prior <- LogisticNormalMixture(
    comp1 = comp1,
    comp2 = comp2,
    weightpar = c(a = 1, b = 1),
    ref_dose = 50
)
normal_mix_prior


## -----------------------------------------------------------------------------
normal_mix_fix_weights_prior <- LogisticNormalFixedMixture(
    components = list(
        comp1 = comp1,
        comp2 = comp2
    ),
    weights = c(0.5, 0.5),
    ref_dose = 50
)
normal_mix_fix_weights_prior


## -----------------------------------------------------------------------------
body(normal_mix_prior@priormodel)
body(normal_mix_fix_weights_prior@priormodel)


## -----------------------------------------------------------------------------
online_mix_prior <- LogisticLogNormalMixture(
    share_weight = 0.5, # probability that the current and external data are from the same distribution
    ref_dose = 50,
    mean = comp1@mean,
    cov = comp1@cov
)


## -----------------------------------------------------------------------------
#| message: false
data_share <- DataMixture(
    doseGrid = c(25, 50, 100, 200, 300),
    x = c(25, 25, 50, 50),
    y = c(0, 1, 0, 1),
    xshare = c(25, 25, 50, 50, 100, 100, 300, 300),
    yshare = c(0, 0, 0, 1, 0, 0, 0, 1)
)


## -----------------------------------------------------------------------------
body(online_mix_prior@priormodel)
body(online_mix_prior@datamodel)


## -----------------------------------------------------------------------------
samples_share <- mcmc(data_share, online_mix_prior, McmcOptions())
plot(samples_share, online_mix_prior, data_share)


## -----------------------------------------------------------------------------
mean(samples_share@data$comp == 2)


## -----------------------------------------------------------------------------
model <- TITELogisticLogNormal(
    mean = c(1.33, 1.49),
    cov = matrix(c(1.826, 0.0209, 0.0209, 0.0245), nrow = 2),
    ref_dose = 56,
    weight_method = "linear"
)
increments <- IncrementsRelative(intervals = c(0, 20), increments = c(2, 1))
next_best <- NextBestMTD(
    # Note that we use a "nonparametric" rule here.
    target = 0.3,
    derive = function(mtd_samples) mean(mtd_samples)
)
stopping <- StoppingMinPatients(n = 20)
size <- CohortSizeConst(3)
t_max <- 42
# Note that we need to use the DataDA class because we need to have
# DLT times in the data.
empty_data <- DataDA(doseGrid = seq(from = 2, to = 50, by = 2), Tmax = t_max)
safety_window <- SafetyWindowConst(gap = c(7, 5), follow = 7, follow_min = 14)


## -----------------------------------------------------------------------------
truth <- function(x) {
    plogis(2 + 3 * log(x / 56))
}

onset <- 3
exp_cond_cdf <- function(x) {
    1 -
        (pexp(x, 1 / onset, lower.tail = FALSE) -
            pexp(t_max, 1 / onset, lower.tail = FALSE)) /
            pexp(t_max, 1 / onset)
}


## -----------------------------------------------------------------------------
curve(exp_cond_cdf, from = 0, to = t_max, ylab = "CDF", xlab = "Time")
abline(v = onset)


## -----------------------------------------------------------------------------
design <- DADesign(
    model = model,
    increments = increments,
    nextBest = next_best,
    stopping = stopping,
    cohort_size = size,
    data = empty_data,
    safetyWindow = safety_window,
    startingDose = 2
)
sims <- simulate(
    object = design,
    truthTox = truth,
    truthSurv = exp_cond_cdf,
    trueTmax = t_max,
    nsim = 100,
    seed = 413,
    mcmcOptions = McmcOptions(
        burnin = 500,
        step = 2,
        samples = 200,
        rng_kind = "Mersenne-Twister",
        rng_seed = 413
    )
)


## -----------------------------------------------------------------------------
summary(sims, truth = truth, target = c(0.29, 0.31))

