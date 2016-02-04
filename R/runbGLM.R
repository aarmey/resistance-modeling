library("MCMCglmm")

data <- read.csv("fullResistanceSet.csv", header = TRUE, sep = ",")

CellLine <- data$CellLine
data$CellLine <- NULL

IDX1 <- CellLine > 2
IDX2 <- CellLine < 3

data1 <- data.frame(data[IDX1,])
CellLine1 <- CellLine[IDX1]

data2 <- data.frame(data[IDX2,])
CellLine2 <- CellLine[IDX2]

N <- 1000

EGFR <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38, data = data1, 
                   nitt = 10*N, thin = N/1000, verbose = FALSE, burnin = N)

HER2 <- MCMCglmm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Akt*Erk, data = data2, 
                   nitt = 10*N, thin = N/1000, verbose = FALSE, burnin = N)

EGFRf <- glm(Viab ~ Akt + Erk + GSK + cJun + JNK + P38 + Akt*Erk, data = data2)

EGFRfnull <- glm(Viab ~ 1, data = data1)

print(EGFRf$aic - EGFRfnull$aic)

rm(data)
rm(data1)
rm(data2)
rm(IDX1)
rm(IDX2)
rm(N)
rm(CellLine)
rm(CellLine1)
rm(CellLine2)