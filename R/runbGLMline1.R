rm(list = ls())

library("MCMCglmm")

dataIn <- read.csv("fullResistanceLine2.csv", header = TRUE, sep = ",")
dataIn$Condition <- NULL
dataIn$CellLine <- NULL

for(i in seq_len(ncol(dataIn))) dataIn[,i] <- scale(dataIn[,i])

N <- 3000

HER2 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn, 
                   nitt = 10*N, verbose = FALSE, burnin = N)

HER2f <- glm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn)
HER2fnull <- glm(Viab ~ 1, data = dataIn)

print(summary(HER2f))
print(HER2f$aic - HER2fnull$aic)

mattx <- as.matrix(HER2$Sol)


panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks = 40)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "black", ...)
}

pairs(mattx, pch='.', diag.panel = panel.hist)