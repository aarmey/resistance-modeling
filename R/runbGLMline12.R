rm(list = ls())

library(MCMCglmm)
library(ggplot2)
library(MASS)
library(coda)
library(emdbook)

files <- c("fullResistanceLine1.csv", "fullResistanceLine2.csv", "fullResistanceLine3.csv", "fullResistanceLineHCC.csv")

loadData <- function (nameIn) {
  dataTemp <- read.csv(nameIn, header = TRUE, sep = ",")
  
  dataTemp$Condition <- NULL
  dataTemp$CellLine <- NULL
  
  for(i in seq_len(ncol(dataTemp))) dataTemp[,i] <- scale(dataTemp[,i])
  
  return(dataTemp)
}

glmmFun <- function (dataIn) {
  N <- 10000
  
  m1 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn, 
                 nitt = 10*N, thin = 10, verbose = FALSE, burnin = N)
  
  fram1 <- as.data.frame(m1$Sol)
  fram1$Model <- rep(1, nrow(fram1))
  
  return (fram1)
}

plotFun <- function (framIn, varsIn, xlimIn, ylimIn) {
  N <- 1000
  
  HPDregionplot(mcmc(framIn[[1]]), vars = varsIn, n = N, prob=0.95, ylim=ylimIn, xlim=xlimIn)
  par(new = TRUE)
  HPDregionplot(mcmc(framIn[[1]]), vars = varsIn, n = N, prob=0.90, ylim=ylimIn, xlim=xlimIn)
  par(new = TRUE)
  HPDregionplot(mcmc(framIn[[2]]), vars = varsIn, n = N, prob=0.95, ylim=ylimIn, xlim=xlimIn)
  par(new = TRUE)
  HPDregionplot(mcmc(framIn[[2]]), vars = varsIn, n = N, prob=0.90, ylim=ylimIn, xlim=xlimIn)
  par(new = TRUE)
  HPDregionplot(mcmc(framIn[[3]]), vars = varsIn, n = N, prob=0.95, ylim=ylimIn, xlim=xlimIn)
  par(new = TRUE)
  HPDregionplot(mcmc(framIn[[3]]), vars = varsIn, n = N, prob=0.90, ylim=ylimIn, xlim=xlimIn)
  par(new = TRUE)
  HPDregionplot(mcmc(framIn[[4]]), vars = varsIn, n = N, prob=0.95, ylim=ylimIn, xlim=xlimIn)
  par(new = TRUE)
  HPDregionplot(mcmc(framIn[[4]]), vars = varsIn, n = N, prob=0.90, ylim=ylimIn, xlim=xlimIn)
}

dataIn <- lapply(files, loadData)
fram <- lapply(dataIn, glmmFun)

plotFun(fram, c("cJun", "Akt"), c(-1, 1), c(-0.5, 2))
