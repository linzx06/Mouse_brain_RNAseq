library(locfdr)
library(edgeR)
source("function_simf_mrf.R")

##########################Simulating the data########################
source("simfunction.R")

####simulation setting 1####
paraMRF=c(-1.1, 2, 0.1, 2, 0.2, 1) ##the mrf parameter
lambda <- 4 ##fold change in the mean
iprob <- 0.5 ##the initial probability of differential expression
dataA <- getsimSeq(option="mrf",paraMRF=paraMRF, b=lambda, iter=3, iprob=iprob)
states <- dataA$states ##the true latent states, 1 means differentially expressed
data <- dataA$data ##read count data
tp <- dataA$tp ##time periods
reg <- dataA$breg ##the layers
sex <- dataA$sex ##sex
group <- dataA$group ##group information for estimating dispersion, samples within the sample group are pooled to estimate dispersion
save(dataA, file="dataA_sim1.rda")

####simulation setting 2####
nDE=0.1 ##the proportion of DE genes
tch=0.5 ##transmission probability in HMM 
perc <- 0.2 ##disturbance from HMM, in the paper, we tried 0.1, 0.2, 0.5
lambda <- 4 ##fold change in the mean
dataA <- getsimSeq(option="hmm", nDE=nDE, tch=tch, perc=perc, b=lambda)
states <- dataA$states ##the true latent states, 1 means differentially expressed
data <- dataA$data ##read count data
tp <- dataA$tp ##time periods
reg <- dataA$breg ##the layers
sex <- dataA$sex ##sex
group <- dataA$group ##group information for estimating dispersion, samples within the sample group are pooled to estimate dispersion
save(dataA, file="dataA_sim2.rda")

#######################Running the algorithm###############
####simulation setting 1####
load("dataA_sim1.rda")
tmp <- calf(data, reg, tp, sex, group, method = "TMM") ##estimate f_0 and f_1
fy1 <- tmp$fy1 ##the non-null density f_1
fy0 <- tmp$fy0 ##the null density f_0
p1 <- tmp$p1 ##probabily of DE, estimated with the R package locfdr
fdrEB <- tmp$fdr ##local fdr, estimated with the R package locfdr, without considering the complex data structure
##run the model
mrfsimf1 <- getMRFDE(fy1=fy1, fy0=fy0, p1=p1, iterEM=200, iterGibbsPost=20000,brPost=10000, skip=5)
mrfsimf1$pfdr1 ##the posterior probability of being differentially expressed

####simulation setting 2
load("dataA_sim2.rda")
tmp <- calf(data, reg, tp, sex, group, method = "TMM") ##estimate f_0 and f_1
fy1 <- tmp$fy1 ##the non-null density f_1
fy0 <- tmp$fy0 ##the null density f_0
p1 <- tmp$p1 ##probabily of DE, estimated with the R package locfdr
fdrEB <- tmp$fdr ##local fdr, estimated with the R package locfdr, without considering the complex data structure
##run the model
mrfsimf2 <- getMRFDE(fy1=fy1, fy0=fy0, p1=p1, iterEM=200, iterGibbsPost=20000,brPost=10000, skip=5)
mrfsimf2$pfdr1 ##the posterior probability of being differentially expressed
