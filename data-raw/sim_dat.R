## Generates the bundled example dataset `sim_dat`.
## A simulated phase-II dose-ranging trial with a binary responder outcome,
## prognostic baseline covariates, and outcomes missing at random.
## Run from the package root:  Rscript data-raw/sim_dat.R
set.seed(2024)

doses <- c(0.00, 0.05, 0.20, 0.40, 0.70, 1.00)   # placebo + 5 active
n_per <- 100                                     # patients per arm
dose  <- rep(doses, each = n_per)
N     <- length(dose)

## prognostic baseline covariates
x1 <- round(rnorm(N), 3)                         # continuous, strongly prognostic
x2 <- rbinom(N, 1, 0.5)                          # binary
x3 <- round(runif(N), 3)                         # continuous on [0, 1]

## binary outcome: prognostic main effect + an Emax dose effect (logit scale)
lp <- -0.4 - 1.2 * x1 + 0.8 * x2 - 1.0 * x3 + 1.3 * dose / (0.2 + dose)
y  <- rbinom(N, 1, plogis(lp))

## outcomes missing at random (~20%), depending on x1, x2
obs <- rbinom(N, 1, plogis(1.5 + 0.6 * x1 - 0.3 * x2))
y[obs == 0] <- NA

sim_dat <- data.frame(dose = dose, y = y, x1 = x1, x2 = x2, x3 = x3)

dir.create("data", showWarnings = FALSE)
save(sim_dat, file = "data/sim_dat.rda", version = 2, compress = "xz")
cat(sprintf("sim_dat: %d rows, %d arms, %.0f%% missing y\n",
            nrow(sim_dat), length(doses), 100 * mean(is.na(sim_dat$y))))
