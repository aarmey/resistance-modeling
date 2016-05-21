rm(list = ls())

library("MCMCglmm")
library("reshape2")

readIn <- function (fileName) {
  dataIn <- read.csv(fileName, header = TRUE, sep = ",")
  dataIn$Condition <- NULL
  dataIn$CellLine <- NULL
  
  for(i in seq_len(ncol(dataIn))) dataIn[,i] <- scale(dataIn[,i])
  
  return(dataIn)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks = 40)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "black", ...)
}

dataIn <- list(readIn("fullResistanceLineHCC.csv"), readIn("fullResistanceLine1.csv"), 
               readIn("fullResistanceLine3.csv"), readIn("fullResistanceLine2.csv"))

CellLines <- c("HCC", "SKBR3", "PC9", "BT474")

N <- 3000

MCMC <- list()
glmModel <- list()
glmNull <- list()

for (i in seq_len(length(dataIn))) {
  MCMC[[i]] <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn[[i]], 
                   nitt = 10*N, verbose = FALSE, burnin = N)
  
  glmModel[[i]] <- glm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = dataIn[[i]])
  glmNull[[i]] <- glm(Viab ~ 1, data = dataIn[[i]])
}

glmModel[[length(dataIn)+1]] <- glm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = rbind(dataIn[[1]], dataIn[[3]]))
glmModel[[length(dataIn)+2]] <- glm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Erk*Akt, data = rbind(dataIn[[2]], dataIn[[4]]))

print(confint.default(glmModel[[5]]))

#print(summary(HER2f))
#print(HER2f$aic - HER2fnull$aic)

#mattx <- as.matrix(HER2$Sol)

#pairs(mattx, pch='.', diag.panel = panel.hist)

trr <- list()

for (i in seq_len(length(dataIn))) {
  names <- attr(MCMC[[i]]$Sol, "dimnames")
  
  trr[[i]] <- data.frame(MCMC[[i]]$Sol)
  
  trr[[i]]["Line"] <- CellLines[[i]]
  trr[[i]]$X.Intercept. <- NULL
  trr[[i]]$JNK <- NULL
  trr[[i]]$P38 <- NULL
  trr[[i]]$Akt.Erk <- NULL
}

matr <- melt(rbind(trr[[1]], trr[[2]], trr[[3]], trr[[4]]),
             id.vars = "Line")


boxplot(value ~ Line * variable, data = matr, 
        col=(c("gold","darkgreen", "grey", "red")), outline = FALSE)
