################################################################################
## Script to prepare data and perform analysis for fruit-feeding butterflies ###
############## data along Trinational Atlantic Forest (TAF) ####################
################################################################################

# 1st: the raw data of occurrence, traits, phylogeny and sampling site information
# were manipulated in the script 1_BasicData 

# 2nd: the data from bioclim and mapbiomas were manipulated in the 2_SpatialData script. 
# In this script we change the pixel resolution (1km2), cut out the variables of 
# interest for the study area, extract the LULC ratio for each sample point and
# for the TAF as a whole.

# 3rd: The occurrence for 258 species at 60 sites was then modeled using the HMSC
# framework, where in addition to the environmental variables (LULC, bioclimate 
# and elevation) we used two functional attributes (FWL and AR), the phylogenetic
# relationships between species and a spatial random effect.

# 4th: we used the parameters estimated by the model to predict the distribution 
# of all species along the TAF, as well as to describe spatial patterns of 
# richness and functional composition.

# 5th: all the images generated in the work (main text and supplementary material)
# are in the cc script.

# Step 1. Setting model structure and fitting the model -------------------

# packages required
library(Hmsc)
library(raster)
library(sp)
library(ape)
library(phytools)

# read data
comm.bfly <- readRDS(here::here("data/processed/comm_objects/comm.bfly.rds"))
dim(comm.bfly)

env.bfly <- readRDS(here::here("data/processed/comm_objects/Envir_covariates.rds"))
dim(env.bfly)

trait.bfly <- readRDS(here::here("data/processed/comm_objects/bfly_trait.rds"))
dim(trait.bfly)

phy.bfly <- readRDS(here::here("data/processed/comm_objects/tree_bfly_AF.rds"))
phy.bfly

# remove the community with NA
which(is.na(env.bfly), arr.ind = T)
comm.bfly <- comm.bfly[-which(rownames(comm.bfly) == "bor1071"),] #bor1002
sort(colSums(comm.bfly))
env.bfly <- env.bfly[-which(rownames(env.bfly) == "bor1071"),] #bor1002
match(rownames(trait.bfly), phy.bfly$tip.label)

# checking the row and col names
match(rownames(comm.bfly), rownames(env.bfly))
match(colnames(comm.bfly), rownames(trait.bfly))

# richness and species prevalence
P <- colMeans(comm.bfly > 0)
S <- rowSums(comm.bfly)

head(sort(P, decreasing = T), n = length(which(P >= 0.7))) # high prevalence (common)
head(sort(P), n = length(which(P <= 0.3))) # low prevalence (rare)

# proportion of species with high, medium and low prevalence
round((length(which(P <= 0.3))*100)/length(P), 2)
round((length(which(P >= 0.7))*100)/length(P), 2)
round((length(which(P > 0.3 & P < 0.7))*100)/length(P), 2)

# FIGURE S1: prevalence and species richness

# correlation between sampling effort and species richness
cor.test(S, env.bfly$Sampling_effort.h.traps.)

# Define the model components
colnames(env.bfly)

# community data
Y <- comm.bfly
head(Y)

# environemntal data
XData <- env.bfly[,c(10:21)]
head(XData)

# define the model structure of fixed effects
XFormula <- ~ Forest_Formation + Open_Formation + 
  Agriculture + Urbanization + Annual_temp + 
  Prec_hot_qtr + elev

# trait data
TrData <- trait.bfly[,c(1,5)]
str(TrData)
sum(is.na(TrData$FWL)); sum(is.na(TrData$AR))

# define the traits formula
TrFormula <- ~FWL + AR

# phylogenetic relationships
phyloTree <- phy.bfly
is.ultrametric(phyloTree)
plot(phyloTree, cex = .5)

# To define a spatial random effect at the level of the site, we need to include the sites ids
# in the studyDesign
studyDesign <- data.frame(Sites = as.factor(env.bfly$Sites_ID))

# define the spatial coordinates for each site
xy <- as.matrix(cbind(long = env.bfly$Longitude, lat = env.bfly$Latitude))
rownames(xy) <- rownames(env.bfly)
head(xy)
rL <- HmscRandomLevel(sData = xy)

# Now we can creates an Hmsc-class object and define the distribution used to model
# the data - a probit model for presence-absence data.
# We define three concurrent models: m.FULL - environmental and spatial components
# m.ENV - only environmental component, and m.SPACE - only spatial component

m.FULL <- Hmsc(Y = Y, XData = XData, XFormula = XFormula, phyloTree = phyloTree, 
               TrData = TrData, TrFormula = TrFormula, distr = "probit",
               studyDesign = studyDesign, ranLevels = list(Sites = rL))
m.FULL

m.ENV <- Hmsc(Y = Y, XData = XData, XFormula = XFormula, phyloTree = phyloTree, 
              TrData = TrData, TrFormula = TrFormula, distr = "probit")
m.ENV

m.SPACE <- Hmsc(Y = Y, XData = XData, XFormula = ~1, phyloTree = phyloTree, 
                TrData = TrData, TrFormula = TrFormula, distr = "probit", 
                studyDesign = studyDesign, ranLevels = list(Sites = rL))
m.SPACE

models <- list(m.FULL, m.ENV, m.SPACE) 

# define some parameters for model run
## for model test
# nChains = 2
# thin = 10
# samples = 100
# transient = 500*thin 
# verbose = 5*thin

# models[[1]] <- sampleMcmc(models[[1]], thin = thin, samples = samples, 
#                             transient = transient, nChains = nChains, verbose = 
#                               verbose)
# saveRDS(models, here::here("output/Model_test.rds"))

nChains = 3
thin = 50
samples = 100
transient = 500*thin 
verbose = 50*thin

for (i in 1:length(models)){
  models[[i]] <- sampleMcmc(models[[i]], thin = thin, samples = samples, 
                            transient = transient, nChains = nChains, verbose = 
                              verbose)
}

saveRDS(models, here::here("output/Model_3ch100samp_21-09.rds"))


# Step 2. Examining MCMC convergence --------------------------------------

# required packages
library(Hmsc)

#m <- readRDS(here::here("output/Model_test.rds"))
m <- readRDS(here::here("output/Model_3ch100samp_21-09.rds"))
m[[1]] # full
m[[2]] # environmental
m[[3]] # spatial

psrf.beta <- psrf.gamma <- psrf.rho <- psrf.alpha <- psrf.omega <- list()

for (i in 1:length(m)) {
  #i=3
  mpost <- convertToCodaObject(m[[i]], spNamesNumbers = c(T,F), 
                                  covNamesNumbers = c(T,F))
  psrf.beta[[i]] <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
  psrf.gamma[[i]] <- gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf
  psrf.rho[[i]] <- gelman.diag(mpost$Rho, multivariate = FALSE)$psrf
  if(i != 2){ # model 2 has only fixed effects, so we cannot access alpha and omega (random effects)
  psrf.alpha[[i]] <- gelman.diag(mpost$Alpha[[1]], multivariate = FALSE)$psrf
  
  # for omega we will take the subset of species pairs
  tmp <- mpost$Omega[[1]] # residual co-occurrence
  sel <- sample(ncol(tmp[[1]]), size=10000)
    for(j in 1:length(tmp)){ 
      tmp[[j]] <- tmp[[j]][,sel]
    }
  psrf.omega[[i]] <- gelman.diag(tmp, multivariate = FALSE)$psrf
  }
}

saveRDS(list(Beta = psrf.beta, Gamma = psrf.gamma, Rho = psrf.rho, Alpha = psrf.alpha, 
                 Omega = psrf.omega), here::here("output/Final_objects/psrf_all.rds"))

# FIG S2: MCMC convergence
# The MCMC convergence diagnostics can be considered satisfactory if for
# most parameters the potential scale reduction factor is close to the ideal value of one.


# Step 3. Evaluating model fit and comparing models -----------------------

m <- readRDS(here::here("output/Model_3ch100samp_21-09.rds"))
# m <- readRDS(here::here("output/Model_test.rds"))

# We start by making a model selection based on WAIC
WAIC.full <- computeWAIC(m[[1]])
WAIC.env <- computeWAIC(m[[2]])
WAIC.space <- computeWAIC(m[[3]])
cbind(full = WAIC.full, env = WAIC.env, space = WAIC.space)

# We next evaluate model fit in terms of explanatory power
# preds have n-species values
preds.full <- computePredictedValues(m[[1]], expected = F)
MF.full <- evaluateModelFit(hM = m[[1]], predY = preds.full)
median(MF.full$RMSE) # accuracy
median(MF.full$AUC) 
median(MF.full$TjurR2) 
sort(MF.full$TjurR2)


preds.env <- computePredictedValues(m[[2]])
MF.env <- evaluateModelFit(hM = m[[2]], predY = preds.env)
median(MF.env$RMSE)
median(MF.env$TjurR2)

preds.space <- computePredictedValues(m[[3]])
MF.space <- evaluateModelFit(hM = m[[3]], predY = preds.space)
median(MF.space$RMSE)
median(MF.space$TjurR2)


# Cross-validation procedure was only employed for the best model

## for model test
# nChains = 2
# thin = 10
# samples = 100
# transient = 500*thin 
# verbose = 5*thin

nChains <- 3
thin <- 50
samples <- 100
transient <- 500*thin 
verbose <- 50*thin
nParallel <- 4
run.cross.validation <- T # start with TRUE when you introduce the script

model.directory <- file.path(".", "output/Final_objects")
filename <- file.path(model.directory, paste0("CVfull_chn_",as.character(nChains),"_samp_",
                                              as.character(samples),"_thin_",as.character(thin)))

if(run.cross.validation){
  partition <- createPartition(m[[1]], nfolds = 4, column = "Sites")
  preds.full <- computePredictedValues(m[[1]], partition = partition, nParallel = nParallel)
  MFCV.full <- evaluateModelFit(hM = m[[1]], predY = preds.full)
  save(partition, MFCV.full, file = filename)
} else {
  load(filename)
}

CV.mean <- c(mean(MF.full$AUC), mean(MFCV.full$AUC), mean(MF.full$TjurR2), 
             mean(MFCV.full$TjurR2))
names(CV.mean) <- c("AUC","AUC (CV)","TjurR2","TjurR2 (CV)")
CV.mean

plot(MF.full$AUC, MFCV.full$AUC)
abline(0,1)

# We observe that the average (over the species) AUC for explanatory power is greater than
# the average AUC for predictive power.
# These results are expected in the sense that explanatory power is typically higher than 
# predictive power, and that AUC and TjurR2 are based on different units, AUC being typically
# higher.


# Step 4. Exploring parameter estimates -----------------------------------

# We proceed with the analyses only for the best model - mFull
# m <- readRDS(here::here("output/Model_3ch100samp_21-09.rds"))
m.FULL <- readRDS(here::here("output/Model_test.rds"))[[1]]

# define some parameters 
nChains <- length(m.FULL$postList)
thin <- m.FULL$thin
samples <- m.FULL$samples
transient <- m.FULL$transient

# We first perform a variance partitioning of the model.
# To be able to group the environmental variables, we look at the design matrix X that Hmsc has constructed
# by applying the XFormula to the XData.

head(m.FULL$X)

# We observe that the columns 2-5 relate to habitat variation and columns 6-7 to climatic variation.
# Arbitrarily, we include the intercept in the habitat variables, and thus perform the variance partitioning
# with the following grouping of the columns of the X matrix.

groupnames <- c("LULC","Bioclim", "elev")
group <- c(1,1,1,1,1,2,2,3)
VP.full <- computeVariancePartitioning(m.FULL, group = group, groupnames = groupnames)
plotVariancePartitioning(m.FULL, VP.full)

# FIG 5: Variance partitioning

# We next construct a beta-plot showing the estimates of species niche parameters
postBeta <- getPostEstimate(m.FULL, parName = "Beta")
plotBeta(m.FULL, post = postBeta, plotTree = F, spNamesNumbers = c(T,F))

# FIG 3: Beta parameters

# In this figure, we observed that species generally responds positivelly to increase in forest coverage
# some species responds positively and another negatively to the annual temperature and elevation 
# 

# We examine next if the species niches are linked to their traits with a Gamma-plot
postGamma <- getPostEstimate(m.FULL, parName = "Gamma")
plotGamma(m.FULL, post = postGamma, supportLevel = 0.9)

# FIG 4: Gamma parameter

# at 0.9 support level traits are linked with species niche

# Another way of examing the influence of traits is to see how much of the variation they
# explain among the responses of the species to their covariates.
VP.full$R2T$Beta
round(sum(VP.full$R2T$Beta), 3)*100
# 56.9% of the variation in species niche is due traits included in model (AR and FWL)

round(VP.full$R2T$Y, 3)*100
# but only 4.2% of the variation in species occurrences are explained by traits

# We next evaluate the posterior distribution of the phylogenetic signal in species niches
mpost.full <- convertToCodaObject(m.FULL)
round(summary(mpost.full$Rho, quantiles = c(0.025, 0.5, 0.975))[[2]],2)
# 2.5%   50% 97.5% 
# 0.33  0.48  0.63

# There is evidence of phylogenetic signal in the species niche since the 
# values did not overlap zero. Furthermore, looking at figure 3 we observe 
# similar responses from phylogenetically close individuals 


# We next illustrate the species associations revealed by the random effects with the corrplot function.

library(corrplot)
OmegaCor <- computeAssociations(m.FULL)
supportLevel <- 0.95
toPlot <- ((OmegaCor[[1]]$support > supportLevel)
           + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *OmegaCor[[1]]$mean
toPlot <- sign(toPlot)
plotOrder <- corrMatOrder(OmegaCor[[1]]$mean, order = "AOE")
corrplot(toPlot[plotOrder, plotOrder], method = "color", tl.cex = 0.5,
         col = colorRampPalette(c("blue", "white", "red"))(255))

# The red and blue colours indicate those species pairs for which the support for
# either a positive or negative association is at least 0.95.
# Note that we have ordered the species (with the function corrMatOrder) so that the cluster of
# co-occurring species is most easily seen. To keep the original species order, write e.g. plotOrder = 1:m$ns

# We next examine at which spatial scale the variation captured by the random effect occurs.
round(summary(mpost.full$Alpha[[1]])$quantiles, 3)
# for Factor 2 there is clear evidence of a spatial signal, as the 95 per cent credible 
# interval does not include zero,with the spatial range of 4km

write.table(round(summary(mpost.full$Alpha[[1]])$quantiles, 3),
            here::here("output/Final_objects/Rnd_effects_alpha.txt"))

# Step 5. Making predictions ----------------------------------------------
# For entire extension of TAF, in a grid cell of 1km2
library(Hmsc)
library(tidyverse)

# m <- readRDS(here::here("output/Model_3ch100samp_21-09.rds"))
m.FULL <- readRDS(here::here("output/Model_3ch100samp_21-09.rds"))[[1]]

# COMMUNITY-LEVEL

# We start by making gradient plots that visualize how the communities vary among the 
# environmental variables, controling the effects of non-focal covariates
head(m.FULL$XData)
m.FULL$XFormula

# Proportion of Forest Formation
Gradient <- constructGradient(m.FULL, focalVariable = "Forest_Formation",
                              non.focalVariables = list("Open_Formation" = list(1), 
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

# Species richness
plotGradient(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)

# traits: index = 2 FWL, index = 3 AR
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 2, showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 3, showData = TRUE)

# Open formations
Gradient <- constructGradient(m.FULL, focalVariable = "Open_Formation",
                              non.focalVariables = list("Forest_Formation" = list(1), 
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

plotGradient(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 2, showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 3, showData = TRUE)

# urbanization
Gradient <- constructGradient(m.FULL, focalVariable = "Urbanization",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

plotGradient(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 2, showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 3, showData = TRUE)

# Agriculture
Gradient <- constructGradient(m.FULL, focalVariable = "Agriculture",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Urbanization" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

plotGradient(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 2, showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 3, showData = TRUE)

# Annual_temp
Gradient <- constructGradient(m.FULL, focalVariable = "Annual_temp",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Prec_hot_qtr" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

plotGradient(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 2, showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 3, showData = TRUE)

# Prec_hot_qtr
Gradient <- constructGradient(m.FULL, focalVariable = "Prec_hot_qtr",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "elev" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

plotGradient(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 2, showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 3, showData = TRUE)

# elevation
Gradient <- constructGradient(m.FULL, focalVariable = "elev",
                              non.focalVariables = list("Forest_Formation" = list(1),
                                                        "Open_Formation" = list(1),
                                                        "Urbanization" = list(1),
                                                        "Agriculture" = list(1),
                                                        "Annual_temp" = list(1),
                                                        "Prec_hot_qtr" = list(1)))
predY <- predict(m.FULL, Gradient = Gradient, expected = TRUE)

plotGradient(m.FULL, Gradient, pred = predY, measure = "S", showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 2, showData = TRUE)
plotGradient(m.FULL, Gradient, pred = predY, measure = "T", index = 3, showData = TRUE)

# SPATIAL PREDICTIONS

# We next perform spatial predictions. See the script 2_M_SpatialData for data manipulation.
# read data - full data

grid <- readRDS(here::here("data/processed/spatial_objects/Clean_data_envir.rds"))
colnames(grid)[3] <- "Annual_temp"
colnames(grid)[4] <- "Prec_hot_qtr"
colnames(grid)[9] <- "Open_Formation"

range(grid$y)
range(grid$x)

# we will create small datasets of the grid to avoid that R "abort" the session

nrow(grid)/10000
1520000 - nrow(grid)
mat.na <- as.data.frame(matrix(NA, nrow = 1520000 - nrow(grid), ncol = ncol(grid)))
colnames(mat.na) <- colnames(grid)

tmp.grid <- rbind(grid, mat.na)
tmp.grid$split <- rep(1:152, each = 10000)
head(tmp.grid)

tmp.grid <- na.omit(tmp.grid)
tmp.grid$split <- as.factor(tmp.grid$split)
split.grid <- split(tmp.grid, tmp.grid$split)

# split.grid have 152 subsets of grids
head(split.grid[[1]])

ggplot(split.grid$`152`, aes(x = x, y = y, fill = Annual_temp)) +
  geom_tile()

# making the prediction for each subset
nParallel <- 3

getS <- function(p){
  return(rowSums(p))
}

for (i in 1:length(split.grid)) {
  #i=1
  tmp <- split.grid[[i]]
  xy.grid <- as.matrix(cbind(tmp$x, tmp$y))
  XData.grid <- as.data.frame(tmp[,match(colnames(m.FULL$XScaled)[-1], colnames(grid))])
  
  Gradient <- prepareGradient(m.FULL, XDataNew = XData.grid, 
                              sDataNew = list(Sites = xy.grid))
  saveRDS(predict(m.FULL, Gradient = Gradient, expected = T, nParallel = nParallel, 
                  predictEtaMean = TRUE), 
          here::here("output", paste("Pred_", i, ".rds", sep = "")))
  saveRDS(Reduce("+", readRDS(here::here("output", paste("Pred_", i, ".rds", sep = ""))))/(300), # 300 is the posterior samples
          here::here("output", paste("Reduced_", i, ".rds", sep = "")))
  saveRDS(simplify2array(lapply(X = readRDS(here::here("output", paste("Pred_", i, ".rds", sep = ""))), FUN = getS)),
          here::here("output", paste("Simplifyed_", i, ".rds", sep = "")))
  
  file <- paste("output/Pred_", i, ".rds", sep = "")
  
  if (file.exists(file)) {
    unlink(file)
    cat("\n This is the subset", i)
  }
}

# read files of mean occupancy probability for each species and sites

grid <- readRDS(here::here("data/processed/spatial_objects/Clean_data_envir.rds"))
colnames(grid)[3] <- "Annual_temp"
colnames(grid)[4] <- "Prec_hot_qtr"
colnames(grid)[9] <- "Open_Formation"
head(grid)


rds_combo <- list.files(path = here::here("output/Reduced"), pattern = "*.rds", full.names = T) %>%
  purrr::map(readRDS) 

for (i in 1:length(rds_combo)) {
  rds_combo[[i]] <- as.data.frame(rds_combo[[i]])
}
memory.size(10^6)
df.mean.occ <- bind_rows(rds_combo)
head(df.mean.occ)

xy <- grid[match(rownames(df.mean.occ), rownames(grid)), 1:2]
plot(xy$x, xy$y)

predY <- as.matrix(df.mean.occ)
colnames(df.mean.occ)

df.mean.occ <- df.mean.occ[,match(m.FULL$phyloTree$tip.label, colnames(df.mean.occ))]

S <- rowSums(predY)
CWM <- (predY %*% m.FULL$Tr)/matrix(rep(S, m.FULL$nt), ncol = m.FULL$nt)

mapData <- data.frame(xy, df.mean.occ, S, CWM, stringsAsFactors=TRUE)
saveRDS(mapData, here::here("data/processed/comm_objects/spatial_pred.rds"))

# Propagating Uncertainty into Predictions

rds_combo_s <- list.files(path = here::here("output/Simplifyed"),
                          pattern = ".rds", full.names = T) %>%
  map(readRDS) 
memory.size(10^6)

for (i in 1:length(rds_combo_s)) {
  rds_combo_s[[i]] <- as.data.frame(rds_combo_s[[i]])
}

df.ric <- bind_rows(rds_combo_s)
head(df.ric)

df.ric$ES <- apply(df.ric, 1, mean) # posterior mean of expected richness
df.ric$sdS <- sqrt(apply(df.ric, 1, var)) # posterior standard deviation

range(df.ric$ES)
range(df.ric$sdS)

df.ric$Latitude <- grid[match(rownames(df.ric), rownames(grid)), "y"]
df.ric$Longitude <- grid[match(rownames(df.ric), rownames(grid)), "x"]

saveRDS(df.ric, here::here("data/processed/comm_objects/Mean_sd_richness.rds"))

# FIG 7 - Spatial predictions
