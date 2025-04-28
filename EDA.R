#install.packages("corrplot")
library(corrplot)
library(splines)

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

plot(abargelie$InfectedBFVivax ~ abargelie$Day, type = "l", ylim = c(0, 200))
lines(abargelie$InfectedBFFalciparum ~ abargelie$Day, col = "blue")
lines(abargelie$InfectedRDTVivax ~ abargelie$Day, col = "green")
lines(abargelie$InfectedRDTFalciparum ~ abargelie$Day, col = "red")

corrFunc = function(test, dat)
{
  woredas = unique(dat$woreda)
  
  corrDF = data.frame(dat[dat$woreda == woredas[1], test])
  colnames(corrDF) = paste0("", woredas[1])
  
  for (w in 2:(length(woredas)))
  {
    tempDF = data.frame(dat[dat$woreda == woredas[w], test])
    colnames(tempDF) = c(paste0("", woredas[w]))
    corrDF = cbind(corrDF, tempDF)
  }
  
  return(corrDF)
}

tempcorr = corrFunc("Blood.film.P..vivax",  dataRaw)
tempcorr = na.omit(tempcorr)
corrDF = data.frame(lapply(tempcorr, as.numeric))
corrMatrix = cor(corrDF)

corrplot(corrMatrix, method = "color")

bsBasisAbergBFVivax = bs(abargelie$Day, df = 20)

fit = lm(abargelie$InfectedBFVivax ~ bsBasisAbergBFVivax)

newx = seq(min(abargelie$Day), max(abargelie$Day), length = 2269)

# Step 4: Predict using new basis
predictedY = predict(fit, newdata = data.frame(x = newx))

# Step 5: Plot
plot(abargelie$Day, abargelie$InfectedBFVivax, pch = 16, main = "B-Spline Fit", xlab = "Day", ylab = "Infected", type = "l")
lines(newx, predictedY, col = "blue", lwd = 2)
