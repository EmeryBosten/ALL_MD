############################################################
# Script that performs Bayesian MD on KU Leuven experiment #
############################################################

############
# set path #
############

setwd("")

#############
# Libraries #
#############

library(readxl)
library(dplyr)
library(tools)
library(rjags)
library(runjags)
library(rstan)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(Metrics)
library(svglite)

###############
## FUNCTIONS ##
###############


# function that fits quadratic model retention model for single compound
# slope = gradient slopes
# rts = retention times of compound
fitQ <- function(slope, rts){
  x = slope
  y = rts
  dataframe = data.frame(y,x)
  
  
  fit <- lm(y ~ x + I(x^2) + I(x^3), data = dataframe)     
  coeff = fit$coefficients
  return(coeff)
}

# function to generate chrom
chrom = function(x,rts,width,Nsubj){
  y = rep(0,length(x))
  for (i in 1:Nsubj){
    y = y + 0.05*dnorm(x,rts[i],width)
  }
  
  return(y)
}

# function to load results of PRIOR
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# function to compute rmse between two vectors
met = function(v1,v2){
  return(rmse(v1,v2))
}

# function to compute resolution
res = function(rt1, rt2){
  return((abs(rt1-rt2)))
}


##########################
## STEP 1 : PRIOR MODEL ##
##########################

##########
## DATA ##
##########

## Load historical data ##
# retention times
rt_struct <- read_excel("./Experimental Setup.xlsx", 
                       sheet = "structuredKULlip")
rt_struct$rt = as.numeric(rt_struct$rt)
rt_struct$S = as.numeric(rt_struct$S)
rt_struct$compound = as.numeric(factor(rt_struct$compound))


# logP and nHBDon
data <- read_excel("./Experimental Setup.xlsx", 
                      sheet = "historical compounds lip")
logP = as.numeric(data$logP)
nHBDon = as.numeric(data$nHBDon)

##############################
## Visualization PRIOR DATA ##
##############################

# distribution retention times 
# Histogram with density plot rt
hrt = ggplot(rt_struct, aes(x=rt)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="blue",binwidth=2.5,alpha=.1)+
  geom_density(alpha=.2, fill="red") 
hrt

ggsave(filename = "rt_dist.svg",path = "./Experiments/Figures/allCompexp/Historical data/", device='svg', dpi=300)

# distribution logP and nHBD
# Histogram with density plot logP
hlogP = ggplot(data, aes(x=logP)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="blue",binwidth=1.5,alpha=.1)+
  geom_density(alpha=.2, fill="red") 
hlogP

ggsave(filename = "logP_dist.svg",path = "./Experiments/Figures/highLogPexp/Historical data/", device='svg', dpi=300)


# barchart nHBDon
hnHBDon = ggplot(data, aes(x=nHBDon)) + 
  geom_bar(colour="black", fill="blue",alpha=.1)
hnHBDon

ggsave(filename = "nHBDon_dist.svg",path = "./Experiments/Figures/allCompexp/Historical data/", device='svg', dpi=300)


# restructure data
rts_gradients<-matrix(rt_struct[order(rt_struct$compound),]$rt,nrow=max(rt_struct$compound),byrow = T)

# compute coefficients
slopes = c(22.5,	18,	15,	12.85,	11.25,	10,	9, 8.18,	7.5, 6.92,	6.43,	6,	5.625,	5.3,	5,	4.7,	4.5)
coeffs = matrix(0, nrow(rts_gradients), 4)
# regression for every compound in list
for (i in 1:nrow(rts_gradients)){
  rts = rts_gradients[i, ]
  coeff = fitQ(slopes, rts)
  coeffs[i, ] = coeff
  }
 
# get mean and standard deviation ~ prior distribution
coefqm = apply(coeffs, 2, mean)
sdqm = apply(coeffs, 2, sd)

## visualization 

r = coefqm[1] + slopes*coefqm[2] + (slopes^2)*coefqm[3] + (slopes^3)*coefqm[4]

plot(rt_struct$S, rt_struct$rt, col = rt_struct$compound, pch = 16)
lines(slopes,r)

# hierarchical model - logP and nHBDon 

coefsA = lm(coeffs[ ,1]~ logP + nHBDon)$coefficients
seA = coef(summary(lm(coeffs[ ,1]~ logP + nHBDon)))[, "Std. Error"]

coefsB = lm(coeffs[ ,2]~ logP + nHBDon)$coefficients
seB = coef(summary(lm(coeffs[ ,2]~ logP + nHBDon)))[, "Std. Error"]

coefsC = lm(coeffs[ ,3]~ logP + nHBDon)$coefficients
seC = coef(summary(lm(coeffs[ ,3]~ logP + nHBDon)))[, "Std. Error"]

coefsD = lm(coeffs[ ,4]~ logP + nHBDon)$coefficients
seD = coef(summary(lm(coeffs[ ,4]~ logP + nHBDon)))[, "Std. Error"]

mlogP = mean(logP)
mnHBDon = mean(nHBDon)

a = coefsA[1] + coefsA[2]*mlogP + coefsA[3]*mnHBDon #sd = 0.61, 0.14, 0.23
b =  coefsB[1] + coefsB[2]*mlogP + coefsB[3]*mnHBDon #sd = 0.038, 0.0088, 0.014
c = coefsC[1] + coefsC[2]*mlogP + coefsC[3]*mnHBDon #sd = 0.00082, 0.00019, 0.0003 
d = coefsD[1] + coefsD[2]*mlogP + coefsD[3]*mnHBDon #sd = 0.00082, 0.00019, 0.0003 

a
b
c
d

coefqm

#Prediction of new compounds 

data <- read_excel("./Experimental Setup.xlsx", 
                   sheet = "Compounds lip")
logP = as.numeric(data$logP)
nHBDon = as.numeric(data$nHBDon)


##############################
## Visualization new data ##
##############################

new_df = data.frame(logP, nHBDon)

# distribution logP and nHBD
# Histogram with density plot logP
hlogP = ggplot(new_df, aes(x=logP)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="blue",binwidth=0.5,alpha=.1)+
  geom_density(alpha=.2, fill="red") 
hlogP

ggsave(filename = "logP_dist.svg",path = "./Experiments/Figures/highLogPexp/New data/", device='svg', dpi=300)


# barchart nHBDon
hnHBDon = ggplot(new_df, aes(x=nHBDon)) + 
  geom_bar(colour="black", fill="blue",alpha=.1)
hnHBDon

ggsave(filename = "nHBDon_dist.svg",path = "./Experiments/Figures/allCompexp/New data/", device='svg', dpi=300)


# Prediction of compound behaviour (retention model) of new compounds based on logP and nHBDon

# predict initial models based on logP and nHBDon values of compounds 
# predicted cubic model parameter values (a, b, c, and d)
a = coefsA[1] + coefsA[2]*logP + coefsA[3]*nHBDon
b = coefsB[1] + coefsB[2]*logP + coefsB[3]*nHBDon
c = coefsC[1] + coefsC[2]*logP + coefsC[3]*nHBDon
d = coefsD[1] + coefsD[2]*logP + coefsD[3]*nHBDon

# range of slopes to predict
predx = seq(4.5, 22.5, by = 0.1)

# predicted initial models
Nsubj = nrow(data)
models = matrix(0,Nsubj,length(predx))
for (i in 1:Nsubj){
  models[i,] = a[i] + b[i]*predx + c[i]*predx^2 + d[i]*predx^3
}

models_gg = data.frame(predx,t(models))

data_long <- melt(models_gg, id = "predx")

g <- ggplot(data_long,            
            aes(x = predx,
                y = value,
                color = variable,
                linetype = variable)) + labs(x = "slope", y = "tr") +  geom_line(linewidth = 0.75) +
  theme(legend.position = "none") 
g

ggsave(filename = "initial_models_new_data.svg",path = "./Experiments/Figures/highLogPexp/New data", device='svg', dpi=300)


######################################
## STEP 2 : Bayesian Model based MD ##
######################################


##############
# STAN MODEL #
##############


# MODEL : relationship between logP and NHBDonor with parameter a, b, and c
# hierarchical model
# rt = a + b*S + c*S^2
# a = a0 + a1*logP + a2*HBDon 
# b = b0 + b1*logP + b2*HBDon 
# c = c0 + c1*logP + c2*HBDon 
# d = d0 + d1*logP + d2*HBDon 
# a0, a1, a2, a, b0, b1, b2, b, c0, c1, c2, d0, d1, and d2 estimated on historical data 
model4String = "
data {
int<lower=0> Ntotal ;
int<lower=1> Nsubj ;
int<lower=1> Npred ;
real x[Ntotal] ;
real y[Ntotal] ;
int s[Ntotal] ;
real p[Nsubj];
real n[Nsubj];
real predx[Npred];
  
// hyperparameters of the prior distributions

real a0mu ;
real a1mu ;
real a2mu ;

real b0mu ;
real b1mu ;
real b2mu ;

real c0mu ;
real c1mu ;
real c2mu ;

real d0mu ;
real d1mu ;
real d2mu ;

real<lower=0> a0sd ;
real<lower=0> a1sd ;
real<lower=0> a2sd ;

real<lower=0> b0sd ;
real<lower=0> b1sd ;
real<lower=0> b2sd ;

real<lower=0> c0sd ;
real<lower=0> c1sd ;
real<lower=0> c2sd ;

real<lower=0> d0sd ;
real<lower=0> d1sd ;
real<lower=0> d2sd ;

real<lower=0> asd ;
real<lower=0> bsd ;
real<lower=0> csd ;
real<lower=0> dsd ;


int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood


}

parameters {

real a0 ;
real a1 ;
real a2 ;

real b0 ;
real b1 ;
real b2 ;

real c0 ;
real c1 ;
real c2 ;

real d0 ;
real d1 ;
real d2 ;

real a[Nsubj] ;
real b[Nsubj] ;
real c[Nsubj] ;
real d[Nsubj] ;

real<lower=0> sigma[Nsubj] ;

}

model {
//initial values are infered from global model 

a0 ~ normal( a0mu  , a0sd ) ; 
a1 ~ normal( a1mu  , a1sd ) ; 
a2 ~ normal( a2mu  , a2sd ) ; 


b0 ~ normal( b0mu  , b0sd ) ; 
b1 ~ normal( b1mu  , b1sd ) ; 
b2 ~ normal( b2mu  , b2sd ) ; 
 

c0 ~ normal( c0mu  , c0sd ) ; 
c1 ~ normal( c1mu  , c1sd ) ; 
c2 ~ normal( c2mu  , c2sd ) ; 

d0 ~ normal( d0mu  , d0sd ) ; 
d1 ~ normal( d1mu  , d1sd ) ; 
d2 ~ normal( d2mu  , d2sd ) ; 

for ( i in 1:Nsubj ) {
a[i] ~ normal( a0 + a1 * p[i] + a2 * n[i] , asd ) ;
b[i] ~ normal( b0 + b1 * p[i] + b2 * n[i] , bsd ) ;
c[i] ~ normal( c0 + c1 * p[i] + c2 * n[i] , csd ) ;
d[i] ~ normal( d0 + d1 * p[i] + d2 * n[i] , dsd ) ;

}
 
sigma ~ normal( 0 , 0.1 ) ;

if(run_estimation == 1){
for ( i in 1:Ntotal ) {
y[i] ~ normal( a[s[i]] + b[s[i]] * x[i] + c[s[i]] * pow(x[i],2) + d[s[i]] * pow(x[i],3) , sigma[s[i]] ) ;
}
}

}

generated quantities {

real y_tilde[Npred,Nsubj]; 

for ( j in 1:Nsubj ) {
  for ( i in 1:Npred ) {
    y_tilde[i,j] = normal_rng( a[j] + b[j] * predx[i] + c[j] * pow(predx[i],2) + d[j] * pow(predx[i],3) , sigma[j]);
}
}

}

" # close quote for modelString

##############
# EXPERIMENT #
##############

# load retention times
# slope, rt, compound name or ID
myData = read_xlsx("./Experimental Setup.xlsx", sheet = "structured KUL lip") # Read the data file; result is a data.frame.
myData$S = as.numeric(myData$S)
myData$compound = as.numeric(factor(myData$compound))
myData$rt = as.numeric(myData$rt)

# reshape format
rts_gradients<-matrix(myData$rt,nrow=max(myData$compound))

## Comparison compound models complete data ##

# generate compound retention models with complete data
# compute model coefficients for every compound
slopes = c(22.5,	18,	15,	12.85,	11.25,	10,	9, 8.18,	7.5, 6.92,	6.43,	6,	5.625,	5.3,	5,	4.7,	4.5)
times = c(4,	5,	6,	7,	8,	9,	10, 11,	12, 13,	14,	15,	16,	17,	18,	19,	20)

coeffs = matrix(0, nrow(rts_gradients), 4)
# regression for every compound in list
for (i in 1:nrow(rts_gradients)){
  rts = rts_gradients[i, ]
  coeff = fitQ(slopes, rts)
  coeffs[i, ] = coeff
}


models_comp = matrix(0,Nsubj,length(predx))
for (i in 1:Nsubj){
  models_comp[i,] = coeffs[i,1] + coeffs[i,2]*predx + coeffs[i,3]*predx^2 + coeffs[i,4]*predx^3
}


# visualize retention models 
models_comp_gg = data.frame(predx,t(models_comp))

data_long <- melt(models_comp_gg, id = "predx")
x = rep(slopes,Nsubj)
y = as.vector(t(rts_gradients[1:Nsubj,]))

color = c(rep("X1",17),rep("X2",17),rep("X3",17),rep("X4",17),rep("X5",17))
preddata = data.frame(y, x, color)
preddata = plyr::rename(preddata, c("y" = "k", "x" = "concentration", "color" = "variable"))
g <- ggplot(data_long,            
             aes(x = predx,
                 y = value,
                 linetype = variable,
                 color = variable)) + labs(x = "slope", y = "tr") +  geom_line(linewidth = 0.75) +
  theme(legend.position = "none") +
  geom_point(data = preddata, aes(x = concentration, y = k, color = variable))
g

ggsave(filename = "models_complete_data.svg",path = "./Experiments/Figures/highLogPexp/New data", device='svg', dpi=300)


## Response surface complete data ##

# compute separation between each compound pair (at each slope data point)
ind = 1
res_mat = matrix(0, 0.5*Nsubj*(Nsubj-1), length(predx))
for (i in 1:(Nsubj-1)){
  for (j in (i+1):Nsubj){
    res_mat[ind,] = res(models_comp[i,],models_comp[j,])
    ind = ind + 1
  }
}

predt = seq(20,4,length.out = length(predx))
res_mat[res_mat>0.15] <- 0.15
# set resolution between 0 and 1
res_mat = res_mat/max(res_mat)
# sum over resolutions for each slope point 
res_tot = (apply(res_mat,2, mean))
res_crit = (apply(res_mat,2, min))
obj_comp = res_tot + res_crit/predt

# plot resolution 
resdata <- data.frame(obj_comp, predx)
resdata <- plyr::rename(resdata, c("obj_comp" = "obj", "predx" = "concentration"))
g <- ggplot(data = resdata, aes(x = concentration, y = obj))  +  geom_line(linetype = "solid", linewidth = 1) +
  labs(x = "slope", y = "obj") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
g

ggsave(filename = "response_surface_all_data.svg",path = "./Experiments/Figures/highLogPexp/New data", device='svg', dpi=300)


##########
# Script #
##########

######################
# DATA PREPROCESSING #
######################

# reshape data into right format for STAN
# gradient times (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
slopes = c(22.5,	18,	15,	12.85,	11.25,	10,	9, 8.18,	7.5, 6.92,	6.43,	6,	5.625,	5.3,	5,	4.7,	4.5)
slope_nums = c(1,17,4) #included experiment
Nsubj = nrow(rts_gradients) # number of compounds
y = c(rts_gradients[1:Nsubj,slope_nums])
x = rep(slopes[slope_nums], each = Nsubj)
s = rep(c(1:Nsubj), length(y)/Nsubj)
logp = logP[1:Nsubj]
nHBD = nHBDon[1:Nsubj]
Ntotal = length(y)

# concentrations to predict
predx = seq(4.5, 22.5, by = 0.1)
Npred = length(predx)
predt = seq(20,4,length.out = Npred)

##priors mean
a0mu = coefsA[1] 
a1mu = coefsA[2]
a2mu = coefsA[3]
b0mu = coefsB[1]
b1mu = coefsB[2]
b2mu = coefsB[3]
c0mu = coefsC[1]
c1mu = coefsC[2]
c2mu = coefsC[3]
d0mu = coefsD[1]
d1mu = coefsD[2]
d2mu = coefsD[3]

# prior sd
a0sd = seA[1]
a1sd = seA[2]
a2sd = seA[3]
b0sd = seB[1]
b1sd = seB[2]
b2sd = seB[3]
c0sd = seC[1]
c1sd = seC[2]
c2sd = seC[3]
d0sd = seD[1]
d1sd = seD[2]
d2sd = seD[3]

asd = sdqm[1]
bsd = sdqm[2]
csd = sdqm[3]
dsd = sdqm[4]

# estimate of likelihood
# to run only prior (without data) = 0
# to run Bayesian inference with data = 1
run_estimation = 1

# MODEL #

dataList = list(
  x = x ,
  y = y ,
  s = s ,
  p = logp ,
  n = nHBD ,
  Ntotal = Ntotal ,
  Nsubj = Nsubj ,
  predx = predx  ,
  Npred = Npred,
  a0mu = a0mu,
  a1mu = a1mu,
  a2mu = a2mu,
  b0mu = b0mu,
  b1mu = b1mu,
  b2mu = b2mu,
  c0mu = c0mu,
  c1mu = c1mu,
  c2mu = c2mu,
  d0mu = d0mu,
  d1mu = d1mu,
  d2mu = d2mu,
  a0sd = a0sd,
  a1sd = a1sd,
  a2sd = a2sd,
  b0sd = b0sd,
  b1sd = b1sd,
  b2sd = b2sd,
  c0sd = c0sd,
  c1sd = c1sd,
  c2sd = c2sd,
  d0sd = d0sd,
  d1sd = d1sd,
  d2sd = d2sd,
  asd = asd, 
  bsd = bsd, 
  csd = csd, 
  dsd = dsd, 
  run_estimation = run_estimation
)

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code=model4String ) 


# RAM to estimate your model in parallel
options(mc.cores = parallel::detectCores())

# automatically save a bare version of a compiled Stan program
rstan_options(auto_write = TRUE)

#Inferred parameters
parameters = c( "a0", "a1", "a2", "a" , "b0", "b1", "b2", "b" , "c0", "c1", "c2", "c", "d0", "d1", "d2", "d", "sigma", "y_tilde")

######################
# Bayesian inference #
######################

# set sampling parameters
adaptSteps = 4000  # Number of steps to "tune" the samplers
burnInSteps = 4000 # Number of steps for Burn-In
nChains = 3 
thinSteps = 1
numSavedSteps=15000

# Get MC sample of posterior:
stanFit <- sampling( object=stanDso , 
                     data = dataList , 
                     pars = parameters ,
                     chains = nChains ,
                     iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                              +burnInSteps ) , 
                     warmup = burnInSteps , 
                     #init = initsList , # optional
                     thin = thinSteps )
# For consistency with JAGS-oriented functions in DBDA2E collection, 
# convert stan format to coda format:
codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                 function(x) { mcmc(as.array(stanFit)[,x,]) } ) )


# save results
save(codaSamples, file = "coda.samples.KUL.prior.lip_cmps.RData")
save(stanFit, file = "stanFit.KUL.prior.lip_cmps.RData")

# Load results
load("coda.samples.KUL.prior.comp_cmps.RData")
load("stanFit.KUL.prior.comp_cmps.RData")

# Load PRIOR results
stanFitPrior <- loadRData("stanFit.KUL.prior.lip_cmps.RData")

############################
# Visualization of results #
############################

# plot which shows density function of predicted (mean) retention times across all compounds
y_tilde <- as.matrix(stanFit, pars = "y_tilde")
mrt = colMeans(y_tilde)
df = data.frame(predRt = mrt)
p <- ggplot(df, aes(x=predRt)) + 
  geom_density(color="darkblue", fill="lightblue") + geom_vline(aes(xintercept=median(predRt)),
                                                                color="blue", linetype="dashed", linewidth=1)
p

source("C:/Users/emery/OneDrive - KU Leuven/Additional reading/Bayesian theory/DBDA2Eprograms/DBDA2Eprograms/DBDA2E-utilities.R")

# Diagnostic plots
# coefficients a
diagMCMC( codaObject=codaSamples , parName="a0", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png"  )
diagMCMC( codaObject=codaSamples , parName="a1", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png")
diagMCMC( codaObject=codaSamples , parName="a2", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
# coefficients b
diagMCMC( codaObject=codaSamples , parName="b0", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
diagMCMC( codaObject=codaSamples , parName="b1", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
diagMCMC( codaObject=codaSamples , parName="b2", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
# coefficients c
diagMCMC( codaObject=codaSamples , parName="c0", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
diagMCMC( codaObject=codaSamples , parName="c1", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
diagMCMC( codaObject=codaSamples , parName="c2", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )

# coefficients d
diagMCMC( codaObject=codaSamples , parName="d0", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
diagMCMC( codaObject=codaSamples , parName="d1", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
diagMCMC( codaObject=codaSamples , parName="d2", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )

# retention model coefficients compound 1
diagMCMC( codaObject=codaSamples , parName="a[1]", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
diagMCMC( codaObject=codaSamples , parName="b[1]", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
diagMCMC( codaObject=codaSamples , parName="c[1]", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )
diagMCMC( codaObject=codaSamples , parName="d[1]", saveName ="./Experiments/Figures/highLogPexp/Experiment 1/", saveType = "png" )


## density plots retention model coefficients a, b, and c for compound 1
#coefficient a
a1 <- as.matrix(stanFit, pars = "a[5]") #posterior distribution for coefficient a of compound 1
a1p = as.matrix(stanFitPrior, pars = "a[5]") #prior distribution for coefficient a of compound 1
df = data.frame(a_cmp1 = a1[,1], a_prior_cmp1 = a1p[,1])
data = melt(df)
# plot density of posterior distribution only
da1 <- ggplot(df, aes(x=a_cmp1)) + 
  geom_density(color="darkblue", fill="lightblue") + geom_vline(aes(xintercept=median(a_cmp1)),
                                                                color="blue", linetype="dashed", linewidth=1)
da1

ggsave(filename = "dist_a_cmp1.svg",path = "./Experiments/Figures/allCompexp/New data/", device='svg', dpi=300)

# plot density of prior and posterior distribution overlayed
ggplot(data, aes(x=value, fill=variable)) +
  geom_density(alpha=.25) + labs(y= "density", x = "a_cmp6")

ggsave(filename = "dist_a_overlay_cmp5.svg",path = "./Experiments/Figures/highLogPexp/Experiment 0/", device='svg', dpi=300)


#coeffficent b
b1 <- as.matrix(stanFit, pars = "b[5]")
b1p = as.matrix(stanFitPrior, pars = "b[5]") #prior distribution for coefficient a of compound 1
df = data.frame(b_cmp1 = b1[,1], b_prior_cmp1 = b1p[,1])
data = melt(df)
# plot density of posterior distribution only
db1 <- ggplot(df, aes(x=b_cmp1)) + 
  geom_density(color="darkblue", fill="lightblue") + geom_vline(aes(xintercept=median(b_cmp1)),
                                                                color="blue", linetype="dashed", linewidth=1)
db1

ggsave(filename = "dist_b_cmp6.svg",path = "./Experiments/Figures/allCompexp/New data/", device='svg', dpi=300)


# plot density of prior and posterior distribution overlayed

ggplot(data, aes(x=value, fill=variable)) +
  geom_density(alpha=.25) + labs(y= "density", x = "b_cmp1")

ggsave(filename = "dist_b_overlay_cmp5.svg",path = "./Experiments/Figures/highLogPexp/Experiment 0/", device='svg', dpi=300)


#coefficient c
c1 <- as.matrix(stanFit, pars = "c[5]")
c1p = as.matrix(stanFitPrior, pars = "c[5]") #prior distribution for coefficient a of compound 1
df = data.frame(c_cmp1 = c1[,1], c_prior_cmp1 = c1p[,1])
data = melt(df)
# plot density of posterior distribution only
dc1 <- ggplot(df, aes(x=c_cmp1)) + 
  geom_density(color="darkblue", fill="lightblue") + geom_vline(aes(xintercept=median(c_cmp1)),
                                                                color="blue", linetype="dashed", linewidth=1)
dc1

ggsave(filename = "dist_c_cmp6.svg",path = "./Experiments/Figures/allCompexp/New data/", device='svg', dpi=300)

# plot density of prior and posterior distribution overlayed

ggplot(data, aes(x=value, fill=variable)) +
  geom_density(alpha=.25) + labs(y= "density", x = "c_cmp1")

ggsave(filename = "dist_c_overlay_cmp5.svg",path = "./Experiments/Figures/highLogPexp/Experiment 0/", device='svg', dpi=300)

#coefficient d
d1 <- as.matrix(stanFit, pars = "d[5]")
d1p = as.matrix(stanFitPrior, pars = "d[5]") #prior distribution for coefficient a of compound 1
df = data.frame(d_cmp1 = d1[,1], d_prior_cmp1 = d1p[,1])
data = melt(df)
# plot density of posterior distribution only
dd1 <- ggplot(df, aes(x=d_cmp1)) + 
  geom_density(color="darkblue", fill="lightblue") + geom_vline(aes(xintercept=median(d_cmp1)),
                                                                color="blue", linetype="dashed", linewidth=1)
dd1

ggsave(filename = "dist_d_cmp6.svg",path = "./Experiments/Figures/allCompexp/New data/", device='svg', dpi=300)

# plot density of prior and posterior distribution overlayed

ggplot(data, aes(x=value, fill=variable)) +
  geom_density(alpha=.25) + labs(y= "density", x = "d_cmp1")

ggsave(filename = "dist_d_overlay_cmp5.svg",path = "./Experiments/Figures/highLogPexp/Experiment 0/", device='svg', dpi=300)

###############
# PREDICTIONS #
###############


## Plot posterior predictive distributions of analyte retention factor

for (i in 1:Nsubj){
  
  n = i # compound
  x = slopes
  y = rts_gradients[n,]
  xdata <- data.frame(y, x)
  xdata <- plyr::rename(xdata, c("y" = "k", "x" = "concentration"))
  
  x = slopes[slope_nums]
  y = rts_gradients[n,slope_nums]
  preddata = data.frame(y, x)
  preddata = plyr::rename(preddata, c("y" = "k", "x" = "concentration"))
  
  pred <- as.data.frame(stanFit, pars = "y_tilde")[,((n-1)*length(predx)+1):(n*length(predx))] %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.05),
              median = quantile(value, probs = 0.5),
              ub = quantile(value, probs = 0.95))
  
  g1 <- ggplot() + geom_point(data = xdata, aes(x = concentration, y = k)) +
    labs(x = "slope", y = "tr") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) + geom_point(data = preddata, aes(x = concentration, y = k), colour = "red") + geom_line(data = pred, aes(x = predx, y = median)) +
    geom_ribbon(aes(predx, ymin = pred$lb, ymax = pred$ub), alpha = 0.25)
  g1
  
  filename = gsub(' ','',paste("nexp3_cmp",as.character(n),".svg"))
  filename
  ggsave(filename = filename,path = "./Experiments/Figures/highLogPexp/Experiment 0/", device='svg', dpi=300)
}

# median parameter values 
a = apply(as.data.frame(stanFit, pars = "a"),2,mean) 
b = apply(as.data.frame(stanFit, pars = "b"),2,mean) 
c = apply(as.data.frame(stanFit, pars = "c"),2,mean)
d = apply(as.data.frame(stanFit, pars = "d"),2,mean) 

# create compound models 
models = matrix(0,Nsubj,length(predx))
for (i in 1:Nsubj){
  # transform back with exponential
  models[i,] = a[i] + b[i]*predx + c[i]*predx^2 + d[i]*predx^3
}

# visualize retention models 
models_gg = data.frame(predx,t(models))

data_long <- melt(models_gg, id = "predx")
x = rep(slopes[slope_nums],Nsubj)
y = as.vector(t(rts_gradients[1:Nsubj,slope_nums]))
preddata = data.frame(y, x)
preddata = plyr::rename(preddata, c("y" = "k", "x" = "concentration"))

g2 <- ggplot(data_long,            
             aes(x = predx,
                 y = value,
                 color = variable)) + labs(x = "slope", y = "tr") +  geom_line(aes(linetype=variable), linewidth = 0.75) +
  theme(legend.position = "none") +
  geom_point(data = preddata, aes(x = concentration, y = k), color = "red")
g2

ggsave(filename = "nexp0_models_combined.svg",path = "./Experiments/Figures/highLogPexp/New data/", device='svg', dpi=300)


# correlation coefficient with complete models
mets = matrix(0,nrow(models),2)
for (i in 1:nrow(models)){
  mets[i,1] = met(models_comp[i,],models[i,])
  mets[i,2] = cor(models_comp[i,],models[i,])^2
}
mets

#####################
## Active Learning ##
#####################

# compute average Change of variance for entire domain of x for each xs 
y_tilde_mat <- as.matrix(stanFit, pars = "y_tilde")

matvarchange = matrix(0,Nsubj,length(predx))

for (n in 0:Nsubj-1) {
  y_tilde = y_tilde_mat[,(n*length(predx)+1):((n+1)*length(predx))]
  vecvarchange = rep(NA, length(predx))
  for (j in 1:length(predx)) {
    xs = y_tilde[,j]
    sumdeltavar = 0
    for (i in 1:length(predx)) {
      xr = y_tilde[,i]
      deltavarxr = cov(xr,xs)^2*(1+var(xs))
      sumdeltavar = sumdeltavar + deltavarxr
    }
    avdeltavar = sumdeltavar/length(predx)
    vecvarchange[j] = avdeltavar
    
  }
  matvarchange[n+1,] = vecvarchange
}

#mean change of variance across all analytes 
meanvarchange = apply(matvarchange,2, mean)

vardata <- data.frame(meanvarchange, predx)
vardata <- plyr::rename(vardata, c("meanvarchange" = "deltaVariance", "predx" = "concentration"))
g3 <- ggplot() + geom_line(data = vardata, aes(x = concentration, y = deltaVariance),linetype = "twodash", linewidth = 1.5) +
  labs(x = "slope", y = "delta var") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8))
g3

which.max(vecvarchange) 

# compute separation between each compound pair (at each slope data point)
ind = 1
res_mat = matrix(0, 0.5*Nsubj*(Nsubj-1), length(predx))
for (i in 1:(Nsubj-1)){
  for (j in (i+1):Nsubj){
    res_mat[ind,] = res(models[i,],models[j,])
    ind = ind + 1
  }
}

# cap on resolution
res_mat[res_mat>0.15] <- 0.15
# set resolution between 0 and 1
res_mat = res_mat/max(res_mat)
# sum over resolutions for each slope point 
res_tot = (apply(res_mat,2, mean))
res_crit = (apply(res_mat,2, min))
obj = res_tot + res_crit/predt
# plot resolution 

# plot resolution 
resdata <- data.frame(obj, predx)
resdata <- plyr::rename(resdata, c("obj" = "obj", "predx" = "concentration"))
g4 <- ggplot(data = resdata, aes(x = concentration, y = obj))  +  geom_line(linetype = "solid", linewidth = 1) +
  labs(x = "slope", y = "obj") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
g4

ggsave(filename = "exp0_response_surface.svg",path = "./Experiments/Figures/highLogPexp/New data/", device='svg', dpi=300)

which.max(obj)
predx[which.max(obj)]


# combine rate of change and critical resolution to single acquisition function
tun = 1/(1+length(slope_nums))**2
acquisition_fun = tun*(meanvarchange)/(max(meanvarchange)) + (obj/max(obj))

vardata <- data.frame(acquisition_fun, predx)
vardata <- plyr::rename(vardata, c("acquisition_fun" = "deltaVariance", "predx" = "concentration"))
g5 <- ggplot(data = vardata, aes(x = concentration, y = deltaVariance))  +  geom_line(linetype = "solid", linewidth = 1) +
  labs(x = "slope", y = "acquisition") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
g5



figure = ggarrange(g2, g1, g4 + theme(axis.text.y=element_blank()), g5 + theme(axis.text.y=element_blank()), 
                   labels = c("A", "B","C","D"),
                   ncol = 2, nrow = 2)
figure
ggsave(filename = "exp3_summary.svg",path = "./Experiments/Figures/highLogPexp/Experiment 0/", device='svg', dpi=300)


which.max(acquisition_fun)
predx[which.max(acquisition_fun)]
