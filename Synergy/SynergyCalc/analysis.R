library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)

dataTT <- as.data.frame(read.csv("../PC9-U0126-SP6-avg.csv", header = TRUE)) %>%
  melt(id.vars = "X") %>%
  transmute(drugTwo = as.numeric(gsub("[X]","",variable)), drugOne = X, viab = value)

dataTT <- mutate(dataTT, drugOneF = as.factor(drugOne), drugTwoF = as.factor(drugTwo))

#gg <- ggplot(dataTT, aes(x = drugOne, y = viab, color = drugTwo)) + geom_point() + 
#  facet_wrap(~ drugTwo)
#plot(gg)

gg <- ggplot(dataTT, aes(x = drugOneF, y = drugTwoF)) + geom_tile(aes(fill=viab)) + 
  scale_fill_gradient2(low = "red", mid = "white", high = "blue")
plot(gg)


thresh <- filter(dataTT, drugOne < 0.11, drugTwo < 0.11)
thresh <- mean(thresh$viab)/2

drugTwoList <- unique(dataTT$drugTwo)


for (ii in 1:length(drugTwoList)) {
  dataCur <- filter(dataTT, drugTwo == drugTwoList[ii])
  
  
  
  model <- approxfun(dataCur$drugOne, y = dataCur$viab - thresh, ties = mean)
  
  root <- uniroot(model, c(0, 20))
  
  print(drugTwoList[ii])
  print(root$root)
  
  x <- seq(0, 20, by = 0.1)
  
  #plot(x, model(x), col = "tomato", type = "l")
  
}

