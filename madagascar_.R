library(splines)
library(ggplot2)
library(dplyr)
library(lubridate)
library(zoo)
library(tidyr)
library(earth)
library(mgcv)

environmentalData = read.csv("madagascar_cleaned_environmental.csv")
clinicalData = read.csv("madagascar_cleaned_clinical.csv")
metaData = read.csv("madagascar_cleaned_metadata.csv")

clinicalData$X__time = as.Date(clinicalData$X__time)
environmentalData$X__time = as.Date(environmentalData$X__time)

start = min(clinicalData$X__time)

clinicalData$timeSinceStart = as.numeric(clinicalData$X__time - start)

head(environmentalData)
head(clinicalData)
head(metaData)

joinedData = clinicalData %>%
  left_join(environmentalData, by = c("X__time", "commFkt"))

plot(clinicalData$X__time[clinicalData$commFkt == "AMBIABE_AMBALAFASINA"], clinicalData$adjustedCases[clinicalData$commFkt == "AMBIABE_AMBALAFASINA"])
length(clinicalData$X__time[clinicalData$commFkt == "AMBIABE_AMBALAFASINA"])
length(clinicalData$adjustedCases[clinicalData$commFkt == "AMBIABE_AMBALAFASINA"])

subDat = clinicalData[clinicalData$commFkt == "AMBIABE_AMBALAFASINA", ]

adjCases = subDat$adjustedCases
time = subDat$timeSinceStart

bsBasis = bs(time, df = 25)




dailyDates = seq(min(subDat$X__time), max(subDat$X__time), by = "day")
dailyNumeric = as.numeric(dailyDates - start)

dailyBasis = predict(bsBasis, newx = dailyNumeric)

# Using Logs 
m1 = lm(log(adjCases+1) ~ bsBasis) 
predCases = cbind(1, dailyBasis) %*% coef(m1)

predictedCases = pmax(exp(predCases) - 1, 0)

dailyC = data.frame(
  date = dailyDates,
  predictedCases = cbind(1, dailyBasis) %*% coef(m1)
)

plot(clinicalData$X__time[clinicalData$commFkt == "AMBIABE_AMBALAFASINA"], log(clinicalData$adjustedCases[clinicalData$commFkt == "AMBIABE_AMBALAFASINA"]+1), type = "l")
lines(dailyC$date, dailyC$predictedCases, col = "red")

# Not using Logs 

m2 = lm(adjCases ~ bsBasis) 
predCases2 = cbind(1, dailyBasis) %*% coef(m2)

predictedCases2 = pmax(exp(predCases2) - 1, 0)

dailyC2 = data.frame(
  date = dailyDates,
  predictedCases = predCases2
)

plot(clinicalData$X__time[clinicalData$commFkt == "AMBIABE_AMBALAFASINA"], clinicalData$adjustedCases[clinicalData$commFkt == "AMBIABE_AMBALAFASINA"], type = "l")
lines(dailyC2$date, dailyC2$predictedCases, col = "red")

# GAM on Subset

subsetJoined = joinedData[joinedData$commFkt == "AMBIABE_AMBALAFASINA", ]
m3 = gam(adjustedCases ~ s(evi) + s(rainfall) + s(temp) + s(ndvi), data = subsetJoined, family = poisson(link = "log"))

summary(m3)


monthly_env = subsetJoined %>%
  select(X__time, evi, rainfall, temp, ndvi) %>%
  mutate(date = as.Date(X__time))

daily_dates = seq(min(monthly_env$date), max(monthly_env$date), by = "day")

daily_env = data.frame(date = daily_dates)

for (var in c("evi", "rainfall", "temp", "ndvi")) {
  daily_env[[var]] <- approx(
    x = monthly_env$date,
    y = monthly_env[[var]],
    xout = daily_dates,
    rule = 2  # constant extrapolation outside bounds
  )$y
}

daily_env$expected_risk = predict(m3, newdata = daily_env, type = "response")

daily_env = daily_env %>%
  mutate(year = year(date), month = month(date))


monthly_cases = subsetJoined %>%
  mutate(year = year(X__time), month = month(X__time)) %>%
  select(year, month, adjustedCases)

daily_env = left_join(daily_env, monthly_cases, by = c("year", "month"))


daily_env = daily_env %>%
  group_by(year, month) %>%
  mutate(
    predicted_total = sum(expected_risk, na.rm = TRUE),
    scale_factor = adjustedCases / predicted_total,
    estimated_cases = expected_risk * scale_factor
  ) %>%
  ungroup()

ggplot() +
  geom_line(data = daily_env, aes(x = date, y = estimated_cases), color = "red", size = 1) +
  geom_point(data = subsetJoined, aes(x = as.Date(X__time), y = adjustedCases), color = "black", size = 2) +
  geom_line(data = subsetJoined, aes(x = as.Date(X__time), y = adjustedCases), color = "black", linetype = "dashed") +
  labs(
    title = "Daily Estimated Malaria Cases vs. Reported Monthly Totals",
    x = "Date",
    y = "Number of Cases"
  ) +
  theme_minimal()


daily_env %>%
  group_by(year, month) %>%
  summarise(
    sum_estimated_cases = sum(estimated_cases, na.rm = TRUE),
    adjusted_monthly = unique(adjustedCases)
  ) %>%
  mutate(
    diff = sum_estimated_cases - adjusted_monthly,
    percent_error = 100 * diff / (adjusted_monthly + 1e-6)  # Avoid division by 0
  ) -> monthly_check

print(monthly_check, n = 66)


# Heatmap Modelling 

library(sf)
library(osmdata)
library(stringr)
library(fuzzyjoin)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(viridis)
library(ggspatial)
library(prettymapr)
library(ggmap)
library(leaflet)


uniqueAreas = unique(environmentalData$commFkt)

joinedDataClean = na.omit(joinedData)

xSeq = expand.grid("evi" = seq(min(joinedDataClean$evi), max(joinedDataClean$evi), length.out = 50), 
                  "rainfall" = seq(min(joinedDataClean$rainfall), max(joinedDataClean$rainfall), length.out = 50), 
                  "temp" = seq(min(joinedDataClean$temp), max(joinedDataClean$temp), length.out = 50), 
                  "ndvi" = seq(min(joinedDataClean$ndvi), max(joinedDataClean$ndvi), length.out = 50))

regionRisks = data.frame(
  area = uniqueAreas,
  risk = NA_real_,  # explicitly numeric
  stringsAsFactors = FALSE
)


for (i in 1:length(uniqueAreas))
{
  tempSubSet = joinedDataClean[joinedDataClean$commFkt == uniqueAreas[i], ]
  
  if (nrow(tempSubSet) >= 10)
  {
    tempMod = gam(adjustedCases ~ s(evi) + s(rainfall) + s(temp) + s(ndvi), data = tempSubSet, family = poisson(link = "log"))
    
    tempPreds = predict(tempMod, newData = xSeq, type = "response")
    
    tempMeanRisk = mean(tempPreds)
    
    regionRisks$risk[i] = tempMeanRisk
    
    
  }
  
  else
  {
    regionRisks$risk[i] = 0
    
  }
  
}

regionRisks <- regionRisks %>%
  separate(area, into = c("commune", "fokontany"), sep = "_", remove = FALSE) %>%
  mutate(fokontany_clean = str_to_upper(str_replace_all(fokontany, "[^A-Z]", "")))

# Define the bounding box for Madagascar
madagascar_bbox <- getbb("Madagascar")

# Query OSM for villages within the bounding box
madagascar_villages <- opq(bbox = madagascar_bbox) %>%
  add_osm_feature(key = "place", value = "village") %>%
  osmdata_sf()

madagascarVillagesPoints = data.frame(madagascar_villages$osm_points$geometry, madagascar_villages$osm_points$name)
colnames(madagascarVillagesPoints) = c("geometry", "village")

madagascarVillagesPoints$village = toupper(madagascarVillagesPoints$village)
madagascarVillagesPointsFiltered <- madagascarVillagesPoints %>%
  filter(village %in% regionRisks$fokontany) %>%
  distinct(village, .keep_all = TRUE)

madagascarVillagesPointsFiltered <- madagascarVillagesPointsFiltered %>%
  left_join(regionRisks, by = c("village" = "fokontany"))

villagesSf = st_as_sf(madagascarVillagesPointsFiltered)

leaflet(villagesSf) %>%
  addTiles() %>%  # Add default OpenStreetMap base tiles
  addCircleMarkers(
    radius = 5,  # Adjust the radius of the circle markers
    color = ~colorNumeric("YlOrRd", domain = villagesSf$risk)(risk),  # Apply color based on 'risk'
    fillOpacity = 0.7,  # Fill opacity for markers
    popup = ~paste("<strong>Village: </strong>", village, 
                   "<br><strong>Risk: </strong>", risk)  # Popup with village name and risk value
  ) %>%
  addLegend(
    pal = colorNumeric("YlOrRd", domain = villagesSf$risk),  # Color palette for risk values
    values = villagesSf$risk,  # Use 'risk' values for the color scale
    title = "Risk",
    opacity = 1  # Legend opacity
  )
