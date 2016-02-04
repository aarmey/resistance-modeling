rm(list = ls())

library("MCMCglmm")

source("multiplot.R")

dataIn1 <- read.csv("fullResistanceLine1.csv", header = TRUE, sep = ",")
dataIn2 <- read.csv("fullResistanceLine2.csv", header = TRUE, sep = ",")
dataIn3 <- read.csv("fullResistanceLine3.csv", header = TRUE, sep = ",")
dataIn4 <- read.csv("fullResistanceLineHCC.csv", header = TRUE, sep = ",")
dataIn1$Condition <- NULL
dataIn1$CellLine <- NULL
dataIn2$Condition <- NULL
dataIn2$CellLine <- NULL
dataIn3$Condition <- NULL
dataIn3$CellLine <- NULL
dataIn4$Condition <- NULL
dataIn4$CellLine <- NULL

for(i in seq_len(ncol(dataIn1))) dataIn1[,i] <- scale(dataIn1[,i])
for(i in seq_len(ncol(dataIn2))) dataIn2[,i] <- scale(dataIn2[,i])
for(i in seq_len(ncol(dataIn3))) dataIn3[,i] <- scale(dataIn3[,i])
for(i in seq_len(ncol(dataIn4))) dataIn4[,i] <- scale(dataIn4[,i])

N <- 10000

m1 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn1, 
                   nitt = 10*N, thin = 10, verbose = FALSE, burnin = N)

m2 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn2, 
               nitt = 10*N, thin = 10, verbose = FALSE, burnin = N)

m3 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn3, 
               nitt = 10*N, thin = 10, verbose = FALSE, burnin = N)

m4 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn4, 
               nitt = 10*N, thin = 10, verbose = FALSE, burnin = N)

m5 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = rbind(dataIn1, dataIn2), 
               nitt = 10*N, thin = 10, verbose = FALSE, burnin = N)

m6 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = rbind(dataIn3, dataIn4), 
               nitt = 10*N, thin = 10, verbose = FALSE, burnin = N)

m7 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = rbind(dataIn1, dataIn2, dataIn3, dataIn4), 
               nitt = 10*N, thin = 10, verbose = FALSE, burnin = N)

#HER2f <- glm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn)
#HER2fnull <- glm(Viab ~ 1, data = dataIn)

fram1 <- as.data.frame(m1$Sol)
fram1$Model <- rep(1, nrow(fram1))
fram2 <- as.data.frame(m2$Sol)
fram2$Model <- rep(2, nrow(fram2))
fram3 <- as.data.frame(m3$Sol)
fram3$Model <- rep(3, nrow(fram3))
fram4 <- as.data.frame(m4$Sol)
fram4$Model <- rep(4, nrow(fram4))
fram5 <- as.data.frame(m5$Sol)
fram5$Model <- rep(5, nrow(fram5))
fram6 <- as.data.frame(m6$Sol)
fram6$Model <- rep(6, nrow(fram6))
fram7 <- as.data.frame(m7$Sol)
fram7$Model <- rep(7, nrow(fram7))

Mfram <- rbind(fram1, fram2, fram3, fram4)

p1 <- qplot(Mfram$Erk, Mfram$Akt, geom = "density2d", color=factor(Mfram$Model)) + coord_fixed(ratio = 1, xlim = c(-1.5, 1.5))
p2 <- qplot(Mfram$cJun, Mfram$Akt, geom = "density2d", color=factor(Mfram$Model)) + coord_fixed(ratio = 1, xlim = c(-1.5, 1.5))


multiplot(p1, p2, cols = 2)


# panel.hist <- function(x, ...)
# {
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(usr[1:2], 0, 1.5) )
#   h <- hist(x, plot = FALSE, breaks = 40)
#   breaks <- h$breaks; nB <- length(breaks)
#   y <- h$counts; y <- y/max(y)
#   rect(breaks[-nB], 0, breaks[-1], y, col = "black", ...)
# }
# 
# pairs(mattx, pch='.', diag.panel = panel.hist)