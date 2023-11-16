# This file is used to fit the continuous-score occupancy model
## CSOM has been shown to outperform traditional occu model in terms of bias and
## precision (unless annotation level is high - >10% data per site)

library(rlist)
#   If the species is present at the site and calling in the file, 
#   the score is drawn from the Normal distribution:
#     Normal(mu[1], sigma[1])
# 
#   If the species is absent at the site, or present and not calling
#   in the file, the score is drawn from the Normal distribution:
#     Normal(mu[2], sigma[2])

require(nimble)

# Set working directory
# If you are using RStudio and this gives you an error, Source the code instead of Running it
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Set up data and results save directories
# data_save_path <- './simulation_study_data' # Where simulated data were saved
results_save_path <- "./results_CSOM"
if(!dir.exists(results_save_path)){dir.create(results_save_path, recursive = TRUE)}
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)

# skip_duplicates==TRUE will skip re-fitting a model for a scenario
# if there is already a results file for that scenario
skip_duplicates <- TRUE

# NIMBLE custom distribution for the Gaussian mixture distribution of classifier scores:
dContinuousScore <- nimbleFunction(
  run = function(x = double(), theta = double(), mu = double(1), sigma = double(1), z = double(), annotation = double(), log = integer()) {
    if(is.na(annotation)) {
      lp <- (1-z) * dnorm(x, mu[1], sigma[1], log=TRUE) +
        z * log( (1-theta) * dnorm(x, mu[1], sigma[1]) + theta * dnorm(x, mu[2], sigma[2]) )
    } else if(annotation == 0) {
      lp <- (1-z) * dnorm(x, mu[1], sigma[1], log=TRUE) +
        z * (log(1-theta) + dnorm(x, mu[1], sigma[1], log=TRUE))
    } else if(annotation == 1) {
      if(z != 1) stop('error with z value')
      lp <- log(theta) + dnorm(x, mu[2], sigma[2], log=TRUE)
    } else stop('unknown value of annotation')
    returnType(double())
    return(lp)
  }
)

rContinuousScore <- nimbleFunction(
  run = function(n = double(), theta = double(), mu = double(1), sigma = double(1), z = double(), annotation = double()) {
    print('should never call rContinuousScore function')
    returnType(double())
    return(1)
  }
)

registerDistributions(list(
  dContinuousScore = list(
    BUGSdist = 'dContinuousScore(theta, mu, sigma, z, annotation)',
    types = c('mu = double(1)', 'sigma = double(1)')
  )
))

# NIMBLE model code and model object:
code <- nimbleCode({
  
  # Priors for top-level variables
  theta ~ dunif(0, 1)
  psi ~ dunif(0, 1)
  for(i in 1:2) { # Priors for the classifier score distribution parameters
    mu[i] ~ dnorm(0, sd = 10000)
    sigma[i] ~ dunif(0, 10000)
  }
  
  # Model for occupancy
  for(i in 1:NSITES) {
    z[i] ~ dbern(psi)
    for(j in 1:NFILES) {
      score[i, j] ~ dContinuousScore(theta = theta, mu = mu[1:2], sigma = sigma[1:2], z = z[i], annotation = annotation[i, j])
    }
  }
  
  # Assumption to make the two classifier distributions identifiable
  constraint ~ dconstraint(mu[2] >= mu[1])
})


################################################################################

# To use this file with your own data, create a list x 
# save it to an .Rdata file (e.g. "my_filename.Rdata"), unindent the code in the 
# "for" loop below and remove the beginning of the for loop

# list x:
#   NSITES: number of sites total
#   NFILES: number of files per site
#   constraint: 1. Needed for NIMBLE model to enforce mu[1] <= mu[2] (see below)
#   theta_init: 0.5. Initial value to be used for theta in the NIMBLE model.
#   psi_init: 0.5. Initial value to be used for psi in NIMBLE estimation.
#   mu_init (0, 1). Initial values to be used for mu_0 and mu_1 in NIMBLE estimation.
#   sigma_init (1, 2). Initial values to be used for sigma_0 and sigma_1 in NIMBLE estimation.
#   true_score: double() containing the true machine learning score for all files
#   annotation_all: vector of annotations for all annotated files. 
#     Dimensions are (NSITES, NFILES). Each file receives one of the following values:
#     1 (true presence), 0 (true absence), or NA (unannotated)
#   z_data: the z_data (known presences) for the sites in the simulation
#     z_data[i] = 1 if a file from site i has been annotated with a
#     true presence. Otherwise, z_data[i] = NA
#   z_init: Initial values to be used for z (true presence) in NIMBLE estimation.
#     z_init[i] should be set to a number (we recommend 1) for all sites
#     that have a value of NAs in z_data[i] (i.e., true occupancy is unknown)
#     Otherwise, z_init[i] = NA

# Load rds data object with scores
x <- readRDS("C:/Users/ajpc1/Desktop/FIT CSOM TENTSMUIR/Data/x_zero.rds")

# Prior to model fitting, explore score distribution
length(which(x$true_score < 0.001)) / length(x$true_score) # half scores below 0.001 

# Model needs logit scores
x$true_score <- logit(x$true_score)

hist(logit(x$true_score[x$true_score > 0]), main = "Histogram of logit score distributions", xlab = "Logit confidence score")

hist(sample(x$true_score[x$true_score > 0.1], 100), add = TRUE, col = "red")


# add initial values
x$theta_init <- 0.5
x$psi_init <- 0.5
x$mu_init <- c(0, 1)
x$sigma_init <- c(1, 2)
x$z_init <- rep(1, x$NSITES)

counter <- 0

scenario_data <- x

print(scenario_data)

# Garbage collect to clean up memory at beginning of each loop
gc()

# Set up the constants, data, and initial value vectors needed by NIMBLE model
constants <- list(
  NSITES = scenario_data$NSITES,
  NFILES = scenario_data$NFILES
)
data <- list(
  z = scenario_data$z_data,
  annotation = scenario_data$annotation_all,
  score = scenario_data$true_score,
  constraint = 1
)
inits <- list(
  theta = scenario_data$theta_init,
  psi = scenario_data$psi_init,
  mu = scenario_data$mu_init,
  sigma = scenario_data$sigma_init,
  z = scenario_data$z_init
)

# Create the NIMBLE model
Rmodel <- nimbleModel(code = code, constants = constants, data = data, inits = inits)
Rmodel$calculate()

# Build the MCMC, compile, run MCMC
conf <- configureMCMC(Rmodel)

## Note: sites i with any annotation[i,j]=1 will *not* have a sampler
## assigned to z[i]. Those sites have occupancy status fixed at z[i]=1.
## Can check this by uncommenting the line below
## conf$printSamplers()

# Compile NIMBLE model into C for speed
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = TRUE)

# Note: the "annotation" vector (from scenario_data$annotation_all) includes NAs.
# These are placeholders to keep the structure of the annotations and the scores the same
# NIMBLE will notify you of these NAs but this is not a cause for concern.
set.seed(0) #Makes MCMC results replicable
continuous_results <- runMCMC(Cmcmc, niter=10000, nburnin=5000, samples = TRUE, summary = TRUE)

hist(x$true_score, main = "Logit score distribution")

results <- list()
results[[1]] <- continuous_results$summary
################################################################################
##                        INTERPRETING OUTPUT             ######################
################################################################################

