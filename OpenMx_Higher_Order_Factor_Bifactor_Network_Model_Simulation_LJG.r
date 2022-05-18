######################################################################
## Project: Higher-Order Factor, Bifactor, Network Model Simulation
## Date: 2021-06-11
## Author: Lennert J. Groot --- University of Amsterdam
######################################################################

# load required R packages
library( "devtools" )
library( "OpenMx" )
library( "psychonetrics" )
library( "qgraph" )
library( "dplyr" )
library( "MASS")
library( "tictoc" )
library( "ggplot2" )
library( "tidyverse" )
library( "corrr" )

# clear work space, not run
# rm(list = ls())

# Load WAIS-IV correlation matrices and helper functions
# load( url ( "https://github.com/KJKan/mcfarland/blob/master/WAIS_Hungary.Rdata?raw=true" ) )
load( url ( "https://github.com/KJKan/mcfarland/blob/master/WAIS_US.Rdata?raw=true" ) )
load( url ( "https://github.com/KJKan/pame_I/blob/main/data/WAIS_Germany.Rdata?raw=true" ) )
source ( file = "https://raw.githubusercontent.com/KJKan/nwsem/master/nwsem_helper_functions.R" )

# WAIS-IV sample Sizes
n_US <- 1800

# n_Hungary <- 1112
n_Germany <- 1425

# The WAIS-IV subtests; observed variables
yvars <- colnames( WAIS_US )

# Number of observed variables
ny <- length( yvars )

# covariance matrix used in psychonetrics models
cov <- ( n_Germany - 1 )/n_Germany*WAIS_Germany

# latent constructs to be measured (etas)
lvars <- c( "P", # Perceptual
            "V", # Verbal
            "W", # Working Memory
            "S" # Speed
)

# Number of latent constructs
ne <- length( lvars )

# theoretical pattern of factor loadings (simple structure)
lambda_g <- matrix( c( #P V W S
  2, 0, 0, 0, # BD
  0, 2, 0, 0, # SI
  0, 0, 2, 0, # DS
  1, 0, 0, 0, # MR
  0, 1, 0, 0, # VC
  0, 0, 1, 0, # AR
  0, 0, 0, 2, # SS
  1, 0, 0, 0, # VP
  0, 1, 0, 0, # IN
  0, 0, 0, 1, # CD
  0, 0, 1, 0, # LN
  0, 0, 1, 0, # FW
  0, 1, 0, 0, # CO
  0, 0, 0, 1, # CA
  1, 0, 0, 0 # PC
),
ncol = ne,
byrow = TRUE )
lambda_b <- cbind( lambda_g, g = 1 )
lambda_b[ ny, ne + 1 ] <- 2

# Extract a network from the US Sample
NWModel_US <- ggm( covs = ( n_US - 1 )/n_US*WAIS_US,
                   omega = "Full",
                   nobs = n_US )

# Prune it
NWModel_US <- NWModel_US %>% prune( alpha = 0.01, recursive = TRUE )

# Aim for further improvement
NWModel_US <- NWModel_US %>% stepup

# Extract the adjacency matrix and use it as confirmatory network in the Germany sample
omega <- 1*( getmatrix( NWModel_US, "omega" ) !=0 )

# ------------------------------ OpenMx models
# --- Higher order g factor model
# The data as OpenMx object
Data <- mxData( WAIS_Germany, type = "cov", numObs = n_Germany )

# Matrix containing the first order factor loadings
Lambda_g <- mxMatrix( name = 'Lambda',
                      type = 'Full',
                      nrow = ny,
                      ncol = ne,
                      free = lambda_g==1,
                      values = lambda_g/2,
                      labels = label( 'lambda', ny, ne ) )

# Matrix containing second order factor loadings
Gamma <- mxMatrix( name = 'Gamma',
                   type = 'Full',
                   nrow = ne,
                   ncol = 1,
                   free = TRUE,
                   values = 1,
                   labels = label( 'gamma', ne, 1 ) )

# Matrix containing the variance of the general factor
Phi <- mxMatrix( name = 'Phi',
                 type = 'Symm',
                 nrow = 1,
                 free = FALSE,
                 values = 1,
                 labels = 'phi' )

# Matrix containing the residual variances of the first order factors
Psi_g <- mxMatrix( name = 'Psi',
                   type = 'Diag',
                   nrow = ne,
                   free = TRUE,
                   values = diag(ne),
                   labels = label ( 'psi', ne ) )

# Matrix containing the residual variances of the observed variables
Theta <- mxMatrix( name = 'Theta',
                   type = 'Diag',
                   nrow = ny,
                   free = TRUE,
                   values = 1,
                   lbound = 0,
                   labels = label ( 'theta', ny ) )

# The factor model implied variance-covariance matrix
ExpCov_g <- mxAlgebra( name = 'ExpCov',
                       expression = Lambda %*% ( Gamma %*% Phi %*% t( Gamma ) + Psi ) %*% t( Lambda ) +
                         Theta )

# ny by ny identity matrix
I <- mxMatrix( name = 'I',
               type = 'Iden',
               nrow = ny )

# Correlation matrix
ExpCor <- mxAlgebra( name = 'ExpCor',
                     expression = solve( sqrt( I * ExpCov ) ) %*% ExpCov %*% t( solve( sqrt( I * ExpCov ) ) ),
                     dimnames = list( yvars, yvars ) )

# Objective
Obj <- mxExpectationNormal( covariance = 'ExpCov',
                            dimnames = yvars )

# The complete model
gModel <-mxModel( name = 'g Model',
                  Data,
                  Lambda_g,
                  Gamma,
                  Phi,
                  Psi_g,
                  Theta,
                  ExpCov_g,
                  I,
                  ExpCor,
                  Obj,
                  mxFitFunctionML() )

# Fit the model
gFit <- mxRun( gModel )

# Result
( gRes <- summary( gFit, verbose = TRUE ) )

# --- Bifactor model
# Matrix containing the first order factor loadings
Lambda_b <- mxMatrix( name = 'Lambda',
                      type = 'Full',
                      nrow = ny,
                      ncol = ne + 1,
                      free = lambda_b > 0,
                      values = lambda_b/2,
                      labels = label( 'lambda', ny, ne + 1 ) )

# Matrix containing the variances of the factors
Psi_b <- mxMatrix( name = 'Psi',
                   type = 'Diag',
                   nrow = ne + 1,
                   free = FALSE,
                   values = diag( ne + 1 ),
                   labels = label ( 'psi', ne + 1 ) )

# The factor model implied variance-covariance matrix
ExpCov_b <- mxAlgebra( name = 'ExpCov',
                       expression = Lambda %*% Psi %*% t( Lambda ) + Theta )

# The complete model
bModel <-mxModel( name = 'Bifactor model',
                  Data,
                  Lambda_b,
                  Psi_b,
                  Theta,
                  ExpCov_b,
                  I,
                  ExpCor,
                  Obj,
                  mxFitFunctionML() )

# Fit the model
bFit <- mxTryHard( bModel )

# Result
( bRes <- summary( bFit, verbose = TRUE ) )

# compare
mxCompare( bFit, gFit )

# --- Network model
# Matrix containing the scaling parameters
Delta <- mxMatrix( name = 'Delta',
                   type = 'Diag',
                   nrow = ny,
                   ncol = ny,
                   free = TRUE,
                   values = 1,
                   labels = label ('delta', ny ) )

# Matrix containing the partial relations (except the diagonal contains zeroes)
Omega <- mxMatrix( name = 'Omega',
                   type = 'Symm',
                   nrow = ny,
                   ncol = ny,
                   free = omega!=0,
                   values = 0,
                   labels = label( 'omega', ny ) )

# Expected partial correlation Matrix
ExpPcor <- mxAlgebra( name = 'ExpPcor',
                      expression = Omega + I )

# Expected Covariance matrix
ExpCovNW <- mxAlgebra( name = 'ExpCov',
                       expression = Delta %*% solve( I - Omega ) %*% t( Delta ) )

# The complete model
NWModel <- mxModel(name = 'Network Model',
                   Data,
                   I,
                   Delta,
                   Omega,
                   ExpCovNW,
                   ExpPcor,
                   ExpCor,
                   Obj,
                   mxFitFunctionML() )

# Fit the model
NWFit <- mxRun( NWModel )

# Result
( NWRes <- summary( NWFit, verbose = TRUE ) )

# ------------------------------ Compare the OpenMx models
mxCompare( NWFit, bFit )
mxCompare( bFit, gFit )
mxCompare( NWFit, gFit )

# ------------------------------ Extract the model implied correlation matrices
# standardized sigmas
st_sigma_g <- matrix( gFit@output$algebras$`g Model.ExpCor`, ny, dimnames = list( yvars, yvars ) )
st_sigma_b <- matrix( bFit@output$algebras$`Bifactor model.ExpCor`, ny, dimnames = list( yvars, yvars ) )
st_sigma_nw <- matrix( NWFit@output$algebras$`Network Model.ExpCor`, ny, dimnames = list( yvars, yvars
) )
st_sigmas <- list( g = st_sigma_g, b = st_sigma_b, nw = st_sigma_nw )
( gcor_out <- st_sigma_g %>% fashion() )

# write.matrix(gcor_out, "gcor_out.csv", sep = ",")
( bcor_out <- st_sigma_b %>% fashion() )

# write.matrix(bcor_out, "bcor_out.csv", sep = ",")
( nwcor_out <- st_sigma_nw %>% fashion() )

# write.matrix(nwcor_out, "nwcor_out.csv", sep = ",")
# ------------------------------ Simulate Data (OR load simdat_nrep_1000.RData file in line 324)
nrep <- 100 # number of replications per condition
ny <- ncol( st_sigmas[[1]] ) # number of variables
n <- 1425 # sample size used ( German WAIS-IV sample size )
reps <- 1:nrep

# Simulation and model fitting
# generate nrep data sets according to the sigmas
set.seed( 03012021 )
dat <- lapply( st_sigmas,
               function( sigma ) replicate( nrep,
                                            mvrnorm( n,
                                                     rep( 0, ny ),
                                                     sigma ) ) )
gdat <- lapply(reps, function (i) {
  dat[["g"]][,,i]
} )
bdat <- lapply(reps, function (i) {
  dat[["b"]][,,i]
} )
nwdat <- lapply(reps, function (i) {
  dat[["nw"]][,,i]
} )
simdat <- list(gdat = gdat, bdat = bdat, nwdat = nwdat)

# Save Data Array as file
# save(simdat, file = "Simdat_nrep_1000.Rdata")

# Store models in a list
models <- list( HF = gModel, BF = bModel, NW = NWModel )

# Function to fit model x on dataset y
fitMxModel <- function( data, MxModel ) {
  MxModel@data <- NULL
  dat <- mxData( data, type = 'raw' )
  Mu <- mxMatrix( name = 'Mu', nrow = 1, ncol = ncol(data), free = TRUE )
  MxModel@expectation$means <- 'Mu'
  MxModel <- mxModel( MxModel, Mu, dat )
  mxTryHard( MxModel, 25 ) }

# Save sets and reps lists
sets <- list("gdat", "bdat", "nwdat")
FittedModels <- list("HF", "BF", "NW")

#### Fitting models: Run this chunk or load output file (see line 363) from repository to save time as finding solutions will take several hours ####
# Fit all Models on all datasets using lapply
tic("Fitting Models")
fit <- lapply( reps, function (i)
  lapply( sets, function ( set )
    lapply( models, fitMxModel, data = simdat[[set]][[i]] ) ) )

# # Name list elements and sub lists
for (i in seq_along(fit)) {
  names( fit [[i]] ) <- paste(sets)
}

# # Function to get results for all reps, true and fitted models.
GetOut <- function (i) {
  lapply( sets, function ( set )
    lapply( FittedModels, function (model)
      summary(fit[[i]][[set]][[model]], refModels = mxRefModels(fit[[i]][[set]][[model]], run = T) ) ) )
}

# # Gather output
Out <- lapply( reps, GetOut)

# # Name list elements and sub lists in Out
for (i in seq_along(Out)) {
  names( Out [[i]] ) <- paste(sets)
  names( Out [[i]][["gdat"]] ) <- paste(FittedModels)
  names( Out [[i]][["bdat"]] ) <- paste(FittedModels)
  names( Out [[i]][["nwdat"]] ) <- paste(FittedModels)
}
toc()

# Save current List of results for 1000 replications as file
# save(Out, file = "Out_nrep_1000.RData")
#### Optional: Load output file for further analysis ####
# load(url("https://github.com/LJGroot/bifactor/raw/main/Out_nrep_1000.RData"))

#### Analysis of fitted models across conditions ####
# Retrieve specific statistic from specific model summary in replication i
GetStat <- function(i, dat, mod, stat, simplify = "array", USE.NAMES = TRUE) {
  Out[[i]][[dat]][[mod]][[stat]]
}

# Get mean value over all replications for specific fit index for each combination of True/Fitted model
GetFIs <- function (fi) {
  lapply( sets, function ( set )
    lapply( FittedModels, function (model)
      sapply( reps, GetStat, dat = set, mod = model, stat = fi) ) )
}
GetMeanFIs <- function (fi) {
  lapply( sets, function ( set )
    lapply( FittedModels, function (model)
      mean( sapply( reps, GetStat, dat = set, mod = model, stat = fi) ) ) )
}

# Fit Indices of interest
FIs <- list("estimatedParameters", "ChiDoF", "Chi", "p", "fit", "CFI", "TLI", "RMSEA", "RMSEAClose")

# Gather output and rename list elements and sub lists
FIout <- lapply(FIs, GetFIs)
meanFIout <- lapply(FIs, GetMeanFIs)
for (i in seq_along(FIout)) {
  names( FIout ) <- paste(FIs)
}
for (i in seq_along(FIout)) {
  names( FIout [[i]] ) <- paste(sets)
  names( FIout [[i]][["gdat"]] ) <- paste(FittedModels)
  names( FIout [[i]][["bdat"]] ) <- paste(FittedModels)
  names( FIout [[i]][["nwdat"]] ) <- paste(FittedModels)
}
for (i in seq_along(meanFIout)) {
  names( meanFIout ) <- paste(FIs)
}
for (i in seq_along(meanFIout)) {
  names( meanFIout [[i]] ) <- paste(sets)
  names( meanFIout [[i]][["gdat"]] ) <- paste(FittedModels)
  names( meanFIout [[i]][["bdat"]] ) <- paste(FittedModels)
  names( meanFIout [[i]][["nwdat"]] ) <- paste(FittedModels)
}

# Function to Retrieve IC which are nested deeper as element in a [2,3] matrix in Out object
GetIC <- function(i, dat,mod, stat = "informationCriteria", row, col) {
  Out[[i]][[dat]][[mod]][[stat]][row,col]
}

# Get AIC/BIC over all reps for all combinations of True/Fitted models
AICpar <- function (fi) {
  lapply( sets, function ( set )
    lapply( FittedModels, function (model)
      sapply( reps, GetIC, dat = set, mod = model, stat = fi, row = 1, col = 2) ) )
}
AICsample <- function (fi) {
  lapply( sets, function ( set )
    lapply( FittedModels, function (model)
      sapply( reps, GetIC, dat = set, mod = model, stat = fi, row = 1, col = 3) ) )
}
BICpar <- function (fi) {
  lapply( sets, function ( set )
    lapply( FittedModels, function (model)
      sapply( reps, GetIC, dat = set, mod = model, stat = fi, row = 2, col = 2) ) )
}
BICsample <- function (fi) {
  lapply( sets, function ( set )
    lapply( FittedModels, function (model)
      sapply( reps, GetIC, dat = set, mod = model, stat = fi, row = 2, col = 3) ) )
}
IC <- list("informationCriteria")
AICparout <- lapply(IC, AICpar)
AICsampleout <- lapply(IC, AICsample)
BICparout <- lapply(IC, BICpar)
BICsampleout <- lapply(IC, BICsample)

# Combine FI's and IC's to one list and rename unlabeled elements and sub lists
FIT <- c(FIout, AICparout, AICsampleout, BICparout, BICsampleout)
ICs <- list("AIC par", "AIC sample", "BIC par", "BIC sample")
names( FIT ) <- append(FIs, ICs)
for (i in seq_along(FIT) ) {
  names( FIT [[i]] ) <- paste(sets)
  names( FIT [[i]][["gdat"]] ) <- paste(FittedModels)
  names( FIT [[i]][["bdat"]] ) <- paste(FittedModels)
  names( FIT [[i]][["nwdat"]] ) <- paste(FittedModels)
}

# Save Output as file
# save(FIT, file = "FIT_nrep_1000.RData")
# Function to Restructure output to different format to compare fit measures as -> [[FI]][[i]][[set]][[model]]
# instead of [[FI]][[set]][[HF]][i]
CompareFI <- function (fi) {
  lapply(reps, function (i) {
    g <- lapply(FittedModels, function (model) FIT[[fi]][["gdat"]][[model]][i])
    b <- lapply(FittedModels, function (model) FIT[[fi]][["bdat"]][[model]][i])
    n <- lapply(FittedModels, function (model) FIT[[fi]][["nwdat"]][[model]][i])
    x <- list(g,b,n)
    names(x) <- paste(sets)
    for (i in seq_along (x) ) names(x[[i]]) <- paste(FittedModels)
    x
  } )
}

# Collect restructured output for relevant FI's
fi2 <- list( "Chi", "CFI", "TLI", "RMSEA", "fit", "AIC par", "AIC sample", "BIC par", "BIC sample" )
FIT2 <- lapply( fi2, CompareFI ); names(FIT2) <- paste(fi2)

# Functions to established favored models for FI's where desirable values are higher or lower.
PrefMax <- function (fi) {
  lapply(sets, function (set) {
    lapply( reps, function (i) {
      which.max( FIT2[[fi]][[i]][[set]] )
    } )
  } )
}
PrefMin <- function (fi) {
  lapply(sets, function (set) {
    lapply( reps, function (i) {
      which.min( FIT2[[fi]][[i]][[set]] )
    } )
  } )
}
nPrefMax <- function (fi) {
  g <- sum(FIT[[fi]][["gdat"]][["HF"]] >= FIT[[fi]][["gdat"]][["BF"]] & FIT[[fi]][["gdat"]][["NW"]])
  b <- sum(FIT[[fi]][["bdat"]][["BF"]] >= FIT[[fi]][["bdat"]][["HF"]] & FIT[[fi]][["bdat"]][["NW"]])
  nw <- sum(FIT[[fi]][["nwdat"]][["NW"]] >= FIT[[fi]][["nwdat"]][["HF"]] & FIT[[fi]][["nwdat"]][["BF"]])
  n <- list(g,b,nw);names(n) <- paste(FittedModels);n
}
nPrefMin <- function (fi) {
  g <- sum(FIT[[fi]][["gdat"]][["HF"]] <= FIT[[fi]][["gdat"]][["BF"]] & FIT[[fi]][["gdat"]][["NW"]])
  b <- sum(FIT[[fi]][["bdat"]][["BF"]] <= FIT[[fi]][["bdat"]][["HF"]] & FIT[[fi]][["bdat"]][["NW"]])
  nw <- sum(FIT[[fi]][["nwdat"]][["NW"]] <= FIT[[fi]][["nwdat"]][["HF"]] & FIT[[fi]][["nwdat"]][["BF"]])
  n <- list(g,b,nw);names(n) <- paste(FittedModels);n
}

# Generate list of preferred Fitted model for all FI's.
FImax <- list( "CFI", "TLI")
FImin <- list("Chi", "fit", "RMSEA", "AIC par", "AIC sample", "BIC par", "BIC sample" )
PrefMaxout <- lapply(FImax, PrefMax);names (PrefMaxout) <- paste(FImax); for (i in seq_along (PrefMaxout)
) names(PrefMaxout[[i]]) <- paste(sets)
PrefMinout <- lapply(FImin, PrefMin);names (PrefMinout) <- paste(FImin); for (i in seq_along (PrefMinout) )
  names(PrefMinout[[i]]) <- paste(sets)
# Count how many times out of nrep a Fitted Model is preferred when it is also the True Model
PrefMod <- c(PrefMaxout, PrefMinout);for (i in seq_along(PrefMod)) names(PrefMod[[i]]) <- paste(sets)
nPrefMaxout <- lapply(FImax, nPrefMax);names (nPrefMaxout) <- paste(FImax)
nPrefMinout <- lapply(FImin, nPrefMin);names (nPrefMinout) <- paste(FImin)
nPrefMod <- c(nPrefMaxout, nPrefMinout)

# When True Model is NW -> What LV model is preferred, HF | BF?
nwPrefMax <- function (fi) {
  sum(FIT[[fi]][["nwdat"]][["HF"]] < FIT[[fi]][["nwdat"]][["BF"]])
}
( nwPrefMaxout <- lapply(FImax, nwPrefMax) )
nwPrefMin <- function (fi) {
  sum(FIT[[fi]][["nwdat"]][["HF"]] > FIT[[fi]][["nwdat"]][["BF"]])
}
( nwPrefMinout <- lapply(FImin, nwPrefMin) )

# End
sessionInfo()
