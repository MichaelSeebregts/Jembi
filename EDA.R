
dataRaw = read.csv("flattened-environment-clinical-join.csv", header = TRUE)

plot(dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$RDT.P..vivax ~ dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$doy, type = "l", ylim = c(5, 150))
lines(dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$Blood.film.P..vivax ~ dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$doy, col = "blue")
lines(dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$Blood.film.P..falciparum ~ dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$doy, col = "red")
lines(dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$RDT.P..falciparum ~ dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$doy, col = "green")
lines(dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$elev ~ dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == 2012, ]$doy, col = "purple")

head(dataRaw)

unique(dataRaw[dataRaw$woreda == "Abargelie", ]$year)

abargelie_VBF = data.frame("Day" = c(), "InfectedBFVivax" = c(),
                           "InfectedBFFalciparum" = c(), 
                           "InfectedRDTVivax" = c(),
                           "InfectedRDTFalciparum" = c(), 
                           "popAtRisk" = c())

for(i in unique(dataRaw[dataRaw$woreda == "Abargelie", ]$year))
{
  tempdat = dataRaw[dataRaw$woreda == "Abargelie" & dataRaw$year == i,]
  scale = i - 2012
  tempdataFrame = data.frame("Day" = (tempdat$doy + scale*365),
                             "InfectedBFVivax" = tempdat$Blood.film.P..vivax,
                             "InfectedBFFalciparum" = tempdat$Blood.film.P..falciparum,
                             "InfectedRDTVivax" = tempdat$RDT.P..vivax,
                             "InfectedRDTFalciparum" = tempdat$RDT.P..falciparum,
                             "popAtRisk" = tempdat$pop_at_risk)
  abargelie_VBF = rbind(abargelie_VBF, tempdataFrame)
}

plot(abargelie_VBF$InfectedBFVivax ~ abargelie_VBF$Day, type = "l", ylim = c(0, 200))
lines(abargelie_VBF$InfectedBFFalciparum ~ abargelie_VBF$Day, col = "blue")
lines(abargelie_VBF$InfectedRDTVivax ~ abargelie_VBF$Day, col = "green")
lines(abargelie_VBF$InfectedRDTFalciparum ~ abargelie_VBF$Day, col = "red")




regress = function(designMat, y, x)
{
  betaHat = solve(t(designMat)%*%designMat)%*%t(designMat)%*%y
  print(dim(designMat))
  print(length(y))
  
  yHatTemp = designMat%*%betaHat
  
  splineDatTemp = data.frame(cbind(yHatTemp, x))
  return(splineDatTemp)
}

Splines_Custom = function(type, knots, data)
{
  data = data[order(data$x), ]
  splineDat = data.frame("y" = c(), "x" = c())
  if (type == "Piecewise Constant")
  {
    for (i in 1:(length(knots) + 1))
    {
      if (i == 1)
      {
        
        tempPoints = data[data$x < knots[i], ]
        x = tempPoints$x
        y = tempPoints$y
        designMat = matrix(rep(1, length(x)), ncol = 1)
        
        splineDatTemp = regress(designMat, y, x)
        splineDat = splineDatTemp
      }
      
      if(1 < i & i < (length(knots) + 1))
      {
        tempPoints = data[knots[i-1] <= data$x & data$x < knots[i], ]
        x = tempPoints$x
        y = tempPoints$y
        designMat = matrix(rep(1, length(x)), ncol = 1)
        splineDatTemp = regress(designMat, y, x)
        
        splineDat = rbind(splineDat, splineDatTemp)
        
      }
      if(i == (length(knots) + 1))
      {
        tempPoints = data[data$x >= knots[i-1], ]
        x = tempPoints$x
        y = tempPoints$y
        designMat = matrix(rep(1, length(x)), ncol = 1)
        splineDatTemp = regress(designMat, y, x)
        
        splineDat = rbind(splineDat, splineDatTemp)
      }
    }
    
  }
  
  if (type == "Piecewise Linear")
  {
    
    for (i in 1:(length(knots) + 1))
    {
      if (i == 1)
      {
        
        tempPoints = data[data$x < knots[i], ]
        x = tempPoints$x
        y = tempPoints$y
        designMat = matrix(c(rep(1, length(x)), x), ncol = 2)
        
        splineDatTemp = regress(designMat, y, x)
        splineDat = splineDatTemp
      }
      
      if(1 < i & i < (length(knots) + 1))
      {
        tempPoints = data[knots[i-1] <= data$x & data$x < knots[i], ]
        x = tempPoints$x
        y = tempPoints$y
        designMat = matrix(c(rep(1, length(x)), x), ncol = 2)
        splineDatTemp = regress(designMat, y, x)
        
        splineDat = rbind(splineDat, splineDatTemp)
        
      }
      if(i == (length(knots) + 1))
      {
        tempPoints = data[data$x >= knots[i-1], ]
        x = tempPoints$x
        y = tempPoints$y
        designMat = matrix(c(rep(1, length(x)), x), ncol = 2)
        splineDatTemp = regress(designMat, y, x)
        
        splineDat = rbind(splineDat, splineDatTemp)
      }
    }
  }
  
  if (type == "Broken Stick")
  {
    x = data$x
    y = data$y
    designMat = matrix(c(rep(1, length(x)), x), ncol = 2)
    for (i in 1:(length(knots)))
    {
      
      designMat = cbind(designMat, c(pmax(x - knots[i], 0)))
      
    }
    
    splineDatTemp = regress(designMat, y, x)
    splineDat = rbind(splineDat, splineDatTemp)
  }
  
  if (type == "Piecewise Cubic Polynomial")
  {
    
    for (i in 1:(length(knots) + 1))
    {
      if (i == 1)
      {
        
        tempPoints = data[data$x < knots[i], ]
        x = tempPoints$x
        y = tempPoints$y
        designMat = matrix(c(rep(1, length(x)), x, x^2, x^3), ncol = 4)
        
        splineDatTemp = regress(designMat, y, x)
        splineDat = splineDatTemp
      }
      
      if(1 < i & i < (length(knots) + 1))
      {
        tempPoints = data[knots[i-1] <= data$x & data$x < knots[i], ]
        x = tempPoints$x
        y = tempPoints$y
        designMat = matrix(c(rep(1, length(x)), x, x^2, x^3), ncol = 4)
        splineDatTemp = regress(designMat, y, x)
        
        splineDat = rbind(splineDat, splineDatTemp)
        
      }
      if(i == (length(knots) + 1))
      {
        tempPoints = data[data$x >= knots[i-1], ]
        x = tempPoints$x
        y = tempPoints$y
        designMat = matrix(c(rep(1, length(x)), x, x^2, x^3), ncol = 4)
        splineDatTemp = regress(designMat, y, x)
        
        splineDat = rbind(splineDat, splineDatTemp)
      }
    }
    
  }
  
  if (type == "Piecewise Cubic Polynomial and Continous at Knots")
  {
    x = data$x
    y = data$y
    designMat = matrix(c(rep(1, length(x)), x, x^2, x^3), ncol = 4)
    for (i in 1:(length(knots)))
    {
      print(head(designMat))
      
      designMat = cbind(designMat, c(pmax(x - knots[i], 0))^3)
      
    }
    
    splineDatTemp = regress(designMat, y, x)
    splineDat = rbind(splineDat, splineDatTemp)
  }
  
  if (type == "Piecewise Cubic Polynomial and Continous at 1st")
  {
    x = data$x
    y = data$y
    designMat = matrix(c(rep(1, length(x)), x, x^2, x^3), ncol = 4)
    for (i in 1:(length(knots)))
    {
      
      designMat = cbind(designMat, c(pmax(x - knots[i], 0)^3))
      designMat = cbind(designMat, c(pmax(x - knots[i], 0)^2))
      
    }
    
    splineDatTemp = regress(designMat, y, x)
    splineDat = rbind(splineDat, splineDatTemp)
  }
  
  if (type == "Piecewise Cubic Polynomial and Continous at 2nd")
  {
    x = data$x
    y = data$y
    designMat = matrix(c(rep(1, length(x)), x, x^2, x^3), ncol = 4)
    for (i in 1:(length(knots)))
    {
      
      designMat = cbind(designMat, c(pmax(x - knots[i], 0)^3))
      designMat = cbind(designMat, c(pmax(x - knots[i], 0)^2)) 
      designMat = cbind(designMat, c(pmax(x - knots[i], 0))) 
      
    }
    
    splineDatTemp = regress(designMat, y, x)
    splineDat = rbind(splineDat, splineDatTemp)
  }
  
  if (type == "Cubic Spline")
  {
    
  }
  
  if (type == "Natural Spline")
  {
    
  }
  
  if (type == "B Spline Basis")
  {
    
  }
  
  return(splineDat)
}

knots = c(500, 1000, 1500, 2000)

dat = data.frame("y" = abargelie_VBF$InfectedBFVivax, "x" = abargelie_VBF$Day)

splines = Splines_Custom("Piecewise Cubic Polynomial", knots, dat)



