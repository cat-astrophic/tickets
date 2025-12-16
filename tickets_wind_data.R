# This script creates the instrument data

# Loading libraries

library(noaastnr)
library(sf)
library(tigris)

# Project directory

direc <- 'D:/tickets/'

# Get wind data

xxx <- get_stations_info()
xxx <- xxx[which(xxx$state == 'VA'),]
xxx$end_year <- substr(xxx$end, 1, 4)
xxx <- xxx[which(xxx$end_year >= 2022),]

weather_data <- as.data.frame(NULL)

for (i in 1:nrow(xxx)) {
  
  s <- paste(xxx$usaf[i], xxx$wban[i], sep = '-')
  print(s)
  tmp <- get_weather_data(s, 2022)
  tmp$latitude <- rep(xxx$latitude[i], nrow(tmp))
  tmp$longitude <- rep(xxx$longitude[i], nrow(tmp))
  weather_data <- rbind(weather_data, tmp)
  
}

# Make station-day level wind data

weather_data <- weather_data[!is.na(weather_data$wind_dir),]

speeds <- c()
dirs <- c()
stations <- c()
days <- c()
lats <- c()
lons <- c()

weather_data$dates <- substr(weather_data$datetime, 1, 10)

for (s in unique(weather_data$stn)) {
  
  print(s)
  tmp <- weather_data[which(weather_data$stn == s),]
  
  for (d in unique(weather_data$dates)) {
    
    tmpd <- tmp[which(tmp$dates == d),]
    speeds <- c(speeds, mean(tmpd$wind_spd, na.rm = TRUE))
    dirs <- c(dirs, mean(tmpd$wind_dir, na.rm = TRUE))
    stations <- c(stations, s)
    days <- c(days, d)
    lats <- c(lats, tmp$latitude[1])
    lons <- c(lons, tmp$longitude[1])
    
  }
  
}

df <- as.data.frame(cbind(stations, days, speeds, dirs, lats, lons))
colnames(df) <- c('Station', 'Date', 'Speed', 'Direction', 'Latitude', 'Longitude')

# Get fips

co <- counties(state = 'VA')

df2 <- st_as_sf(df, coords = c('Longitude', 'Latitude'))
df2 <- st_set_crs(df2, st_crs(co))

inside <- st_within(df2, co)

hosts <- c()

for (i in 1:nrow(df2)) {
  
  hosts <- c(hosts, co$GEOID[inside[[i]][1]])
  
}

df2$FIPS <- hosts

# Turn this into a regular data frame for saving

dfx <- as.data.frame(df2)
dfx <- dfx[,c(1:4,6)]

# This can now be linked to pollution, cops, and other noaa data

write.csv(dfx, paste0(direc, 'data/wind_data.csv'), row.names = FALSE)

