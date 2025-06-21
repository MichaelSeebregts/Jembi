#install.packages("corrplot")
library(corrplot)
library(splines)

library(waveslim)
library(wavethresh)
library(wavelets)

dataRaw = read.csv("flattened-environment-clinical-join.csv", header = TRUE)

infectedFunc = function(woreda)
{
  dat = data.frame("Day" = c(), "InfectedBFVivax" = c(),
                   "InfectedBFFalciparum" = c(), 
                   "InfectedRDTVivax" = c(),
                   "InfectedRDTFalciparum" = c(), 
                   "popAtRisk" = c())
  
  for(i in unique(dataRaw[dataRaw$woreda == woreda, ]$year))
  {
    tempdat = dataRaw[dataRaw$woreda == woreda & dataRaw$year == i,]
    scale = i - 2012
    tempdataFrame = data.frame("Day" = (tempdat$doy + scale*365),
                               "InfectedBFVivax" = tempdat$Blood.film.P..vivax,
                               "InfectedBFFalciparum" = tempdat$Blood.film.P..falciparum,
                               "InfectedRDTVivax" = tempdat$RDT.P..vivax,
                               "InfectedRDTFalciparum" = tempdat$RDT.P..falciparum,
                               "popAtRisk" = tempdat$pop_at_risk)
    dat = rbind(dat, tempdataFrame)
  }
  return(dat)
}

abargelie = infectedFunc("Abargelie")

plot(abargelie$InfectedBFVivax ~ abargelie$Day, type = "l", ylim = c(0, 200))
lines(abargelie$InfectedBFFalciparum ~ abargelie$Day, col = "blue")
lines(abargelie$InfectedRDTVivax ~ abargelie$Day, col = "green")
lines(abargelie$InfectedRDTFalciparum ~ abargelie$Day, col = "red")


bsBasis = bs(abargelie$Day, df = 30, degree = 3, intercept = FALSE)

m1 = lm(abargelie$InfectedBFVivax ~ bsBasis)
m2 = lm(abargelie$InfectedBFFalciparum ~ bsBasis)
m3 = lm(abargelie$InfectedRDTVivax ~ bsBasis)
m4 = lm(abargelie$InfectedRDTFalciparum ~ bsBasis)




predictedM1 = predict(m1)
predictedM2 = predict(m2)
predictedM3 = predict(m3)
predictedM4 = predict(m4)


plot(abargelie$Day, abargelie$InfectedBFVivax, type = "l")
lines(abargelie$Day, predictedM1, col = "purple")

plot(abargelie$Day, abargelie$InfectedBFFalciparum, type = "l")
lines(abargelie$Day, predictedM2, col = "purple")

plot(abargelie$Day, abargelie$InfectedRDTVivax, type = "l")
lines(abargelie$Day, predictedM3, col = "purple")

plot(abargelie$Day, abargelie$InfectedRDTFalciparum, type = "l")
lines(abargelie$Day, predictedM4, col = "purple")


waveletBFVivax = dwt(abargelie$InfectedBFVivax, n.levels = 4, filter = "la8", boundary = "reflection")
waveletBFFalciparum = dwt(abargelie$InfectedBFFalciparum, n.levels = 4, filter = "la8", boundary = "reflection")
waveletRDTVivax = dwt(abargelie$InfectedRDTVivax, n.levels = 4, filter = "la8", boundary = "reflection")
waveletRDTFalciparum = dwt(abargelie$InfectedRDTFalciparum, n.levels = 4, filter = "la8", boundary = "reflection")



threshold = function(wt, thresh_val)
{
  wt@W <- lapply(wt@W, function(coef) ifelse(abs(coef) > thresh_val, coef, 0))
  return(wt)
}

waveletBFVivaxThresh = threshold(waveletBFVivax, thresh_val = 10)
waveletBFFalciparum = threshold(waveletBFFalciparum, thresh_val = 10)
waveletRDTVivax = threshold(waveletRDTVivax, thresh_val = 10)
waveletRDTFalciparum = threshold(waveletRDTFalciparum, thresh_val = 10)

waveletBFVivaxThreshPred = idwt(waveletBFVivaxThresh)
waveletBFFalciparumThreshPred = idwt(waveletBFFalciparum)
waveletRDTVivaxThreshPred = idwt(waveletRDTVivax)
waveletRDTFalciparumThreshPred = idwt(waveletRDTFalciparum)


plot(abargelie$Day, abargelie$InfectedBFVivax, type = "l")
lines(abargelie$Day, waveletBFVivaxThreshPred, col = "green")

plot(abargelie$Day, abargelie$InfectedBFFalciparum, type = "l")
lines(abargelie$Day, waveletBFFalciparumThreshPred, col = "green")

plot(abargelie$Day, abargelie$InfectedRDTVivax, type = "l")
lines(abargelie$Day, waveletRDTVivaxThreshPred, col = "green")

plot(abargelie$Day, abargelie$InfectedRDTFalciparum, type = "l")
lines(abargelie$Day, waveletRDTFalciparumThreshPred, col = "green")




# corrFunc = function(test, dat)
# {
#   woredas = unique(dat$woreda)
#   
#   corrDF = data.frame(dat[dat$woreda == woredas[1], test])
#   colnames(corrDF) = paste0("", woredas[1])
#   
#   for (w in 2:(length(woredas)))
#   {
#     tempDF = data.frame(dat[dat$woreda == woredas[w], test])
#     colnames(tempDF) = c(paste0("", woredas[w]))
#     corrDF = cbind(corrDF, tempDF)
#   }
#   
#   return(corrDF)
# }
# 
# tempcorr = corrFunc("Blood.film.P..vivax",  dataRaw)
# tempcorr = na.omit(tempcorr)
# corrDF = data.frame(lapply(tempcorr, as.numeric))
# corrMatrix = cor(corrDF)
# 
# corrplot(corrMatrix, method = "color")


