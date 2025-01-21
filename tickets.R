# This script performs the econometric analysis for pollution and tickets

# Loading libraries

library(dplyr)
library(lmtest)
library(sandwich)
library(leaflet)
library(stargazer)
library(tigris)
library(ggplot2)
library(sf)
library(strucchange)
library(MSwM)
library(AER)
library(modelsummary)

# Project directory

direc <- 'D:/tickets/'

# Read in the data

data <- read.csv(paste0(direc, 'data/data.csv'))
counts <- read.csv(paste0(direc, 'data/counts_data.csv'))
crashes <- read.csv(paste0(direc, 'data/CrashData_test_1840667065664721760.csv'))

# Filtering out ambiguous data

data <- data %>% filter(VEHICLE.SEARCHED %in% c('Y', 'N'))
data <- data %>% filter(PERSON.SEARCHED %in% c('Y', 'N'))
data <- data %>% filter(FORCE.USED.BY.OFFICER %in% c('Y', 'N'))
data <- data %>% filter(ENGLISH.SPEAKING %in% c('Y', 'N'))
data <- data %>% filter(GENDER %in% c('FEMALE', 'MALE'))
data <- data %>% filter(ETHNICITY %in% c('NOT HISPANIC OR LATINO', 'HISPANIC OR LATINO'))
data <- data %>% filter(RACE %in% c('WHITE', 'BLACK OR AFRICAN AMERICAN', 'ASIAN OR NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 'AMERICAN INDIAN OR ALASKA NATIVE'))
data <- data %>% filter(RESIDENCY %in% c('RESIDENT OF CITY/COUNTY OF STOP', 'OUT OF STATE RESIDENT', 'OTHER VIRGINIA JURISDICTION RESIDENT'))

# Creating additional variables

data$VEHICLE_SEARCHED <- ifelse(data$VEHICLE.SEARCHED == 'Y', 1, 0)
data$PERSON_SEARCHED <- ifelse(data$PERSON.SEARCHED == 'Y', 1, 0)
data$SEARCHED <- pmax(data$VEHICLE_SEARCHED, data$PERSON_SEARCHED)
data$SEARCHED_BOTH <- (data$VEHICLE_SEARCHED == data$PERSON_SEARCHED) * data$SEARCHED
data$ARRESTED <- ifelse(data$ACTION.TAKEN == 'ARREST', 1, 0)
data$TICKETED <- ifelse(data$ACTION.TAKEN == 'CITATION/SUMMONS', 1, 0)
data$FORCE <- ifelse(data$FORCE.USED.BY.OFFICER == 'Y', 1, 0)

# Create a week of year variable

shit <- function(day, month) {
  
  ref <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  a <- gregexpr('/', day)[[1]][1]
  b <- gregexpr('/', day)[[1]][2]
  
  val <- as.integer(substr(day, a+1, b-1)) + sum(ref[1:(as.integer(month))]) - 31
  
  return(ceiling(val/7))
  
}

weak <- c()

for (i in 1:nrow(data)) {
  
  print(i)
  weak <- c(weak, shit(data$STOP_DATE[i], data$MONTH[i]))
  
}

data$WEEK <- weak

# Create a day of year variable

fuck <- function(day, month) {
  
  ref <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  a <- gregexpr('/', day)[[1]][1]
  b <- gregexpr('/', day)[[1]][2]
  
  val <- as.integer(substr(day, a+1, b-1)) + sum(ref[1:(as.integer(month))]) - 31
  
  return(val)
  
}

day.of.year <- c()

for (i in 1:nrow(data)) {
  
  print(i)
  day.of.year <- c(day.of.year, fuck(data$STOP_DATE[i], data$MONTH[i]))
  
}

data$DAY <- day.of.year

# Compute distance to the mountains (peak of at least 1000 meters) for the IV

options(tigris_use_cache = TRUE)

va <- counties(state = 'VA')
va$GEOID <- as.integer(va$GEOID)

peaks <- read_sf(paste0(direc, '/data/USA_Mountain_Peaks/v107/gsummit.gdb'))
peaks <- st_transform(peaks, 4269)

box <- st_contains(va, peaks)
in.box <- c()

for (i in 1:nrow(box)) {
  
  if (length(box[[i]]) > 0) {
    
    for (j in box[[i]]) {
      
      in.box <- c(in.box, j)
      
    }
    
  }
  
}

peaks <- peaks[in.box,]
peaks <- peaks %>% filter(ELEV_METER >= 1000)

nearest.peak <- st_nearest_feature(va, peaks)
np <- peaks[nearest.peak,]
va$PD <- as.numeric(st_distance(va, np, by_element = TRUE)) / 1000

# Looking at mean daily PM across VA in 2023

day <- c()
dpm <- c()

for (i in 1:max(data$DAY)) {
  
  tmp <- data %>% filter(DAY == i)
  day <- c(day, i)
  dpm <- c(dpm, mean(tmp$PM, na.rm = TRUE))
  
}

fdf <- as.data.frame(cbind(day, dpm))
colnames(fdf) <- c('Day', 'Mean Particulate Matter')

mark <- lm(`Mean Particulate Matter` ~ 1, data = fdf)

msm <- msmFit(mark, k = 2, sw = rep(TRUE, 2))

plotProb(msm, which = 1)

plotProb(msm, which = 2)

plotDiag(msm, regime = 2, which = 3)

fdf$Color <- c(rep('Baseline Pollution', 135), rep('Canadian Wildfires', 70), rep('Baseline Pollution', 160))

ggplot(data = fdf[!is.na(fdf$`Mean Particulate Matter`),], aes(x = Day, y = `Mean Particulate Matter`)) + 
  theme_bw() + 
  ggtitle('Daily Mean Particulate Matter Levels (PM2.5)') + 
  ylab('Micrograms per Cubic Meter') +
  xlab('Day') +
  geom_path(aes(x = Day, y = `Mean Particulate Matter`, group = 1, color = Color), size = 1, alpha = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = 136) +
  geom_vline(xintercept = 206) +
  ylim(c(0,100)) +
  scale_x_continuous(breaks = seq(0, 365, 30), labels = seq(0, 365, 30)) +
  scale_color_manual(breaks = c('Baseline Pollution', 'Canadian Wildfires'), values = c('orange2', 'red4'), name = 'Regime') +
  annotate('text', x = 171, y = 95, label = 'Canadian Wildfires') +
  annotate('text', x = 67, y = 30, label = 'Baseline Pollution') +
  annotate('text', x = 287, y = 30, label = 'Baseline Pollution')

# Compare 2023 to 2022

dpm2 <- c()

for (i in 1:max(data$DAY)) {
  
  tmp <- data %>% filter(DAY == i)
  dpm2 <- c(dpm2, mean(tmp$PM_PY, na.rm = TRUE))
  
}

fdf$PYPM <- dpm2

mark2 <- lm(PYPM ~ 1, data = fdf)

msm2 <- msmFit(mark2, k = 2, sw = rep(TRUE, 2))

plotProb(msm2, which = 1)

plotProb(msm2, which = 2)

plotDiag(msm2, regime = 2, which = 3)

ggplot(data = fdf[!is.na(fdf$`Mean Particulate Matter`),], aes(x = Day, y = `Mean Particulate Matter`)) + 
  theme_bw() + 
  ggtitle('Daily Mean Particulate Matter Levels') + 
  ylab('Micrograms per Cubic Meter') +
  xlab('Day') +
  geom_path(aes(x = Day, y = PYPM, group = 1, color = 'Previous Year Pollution'), size = 2, alpha = 1, color = 'gray') + 
  geom_path(aes(x = Day, y = `Mean Particulate Matter`, group = 1, color = Color), size = 1, alpha = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = 136) +
  geom_vline(xintercept = 205) +
  ylim(c(0,100)) +
  scale_x_continuous(breaks = seq(0, 365, 30), labels = seq(0, 365, 30)) +
  scale_color_manual(breaks = c('Baseline Pollution', 'Canadian Wildfires', 'Previous Year Pollution'), values = c('orange3', 'red4', 'gray'), name = 'Regime') +
  annotate('text', x = 171, y = 95, label = 'Canadian Wildfires') +
  annotate('text', x = 67, y = 30, label = 'Baseline Pollution') +
  annotate('text', x = 287, y = 30, label = 'Baseline Pollution')

# Creating an additional IV with previous year pollution and distance to the mountains

pd <- c()

for (i in 1:nrow(data)) {
  
  print(i)
  pd <- c(pd, va$PD[which(va$GEOID == data$FIPS[i])])
  
}

data$PD <- pd

data$IV_PM_ <- as.integer(data$DAY %in% c(136:205)) * log(data$PD+1)
data$IV_PM10_ <- as.integer(data$DAY %in% c(136:205)) * log(data$PD+1)

data$IV_PM <- (as.integer((data$DAY %in% c(136:205))) * log(data$PD+1)) + data$PM_PY
data$IV_PM10 <- (as.integer((data$DAY %in% c(136:205))) * log(data$PD+1)) + data$PM_PY

# Flagging minors in the data

data$MINOR <- as.integer(data$AGE < 18)

# Prepping for counts regressions

# Updating counts$Date to match data$STOP_DATE

date_fixer <- function (d) {
  
  a <- gregexpr(pattern = '-', d)[[1]][1]
  b <- gregexpr(pattern = '-', d)[[1]][2]
  
  mon <- ifelse(substr(d, 6, 6) == 0, substr(d, 7, 7), substr(d, 6, 7))
  da <- ifelse(substr(d, 9, 9) == 0, substr(d, 10, 10), substr(d, 9, 10))
  
  xxx <- paste0(mon, '/', da, '/', substr(d, 1, 4))
  
  return (xxx)
  
}

counts$DATE2 <- date_fixer(counts$Date)

# Joining

counts <- left_join(counts, data[,which(colnames(data) %in% c('STOP_DATE', 'JURISDICTION', 'PM', 'IV_PM', 'WEEK'))], by = c('DATE2' = 'STOP_DATE', 'Jurisdiction' = 'JURISDICTION'))

counts <- counts[!duplicated(counts), ]

# Combined searches data

counts$Searched <- counts$Vehicles_Searched + counts$Persons_Searched

# Adding counts data to data

data <- left_join(data, counts, by = c('STOP_DATE' = 'DATE2', 'AGENCY.NAME' = 'Agency', 'JURISDICTION' = 'Jurisdiction'))

# Update some column names

colnames(data)[which(colnames(data) == 'PM.x')] <- 'PM'
colnames(data)[which(colnames(data) == 'IV_PM.x')] <- 'IV_PM'
colnames(data)[which(colnames(data) == 'WEEK.x')] <- 'WEEK'

# Adding another IV

data$XXX <- data$PD * as.integer(data$DAY %in% c(136:205))

# Adding a day-of-week variable

data$DOW <- data$DAY %% 7

# Previous day pollution data

data$DOY <- (data$WEEK - 1)*7 + data$DOW
pd.pm <- c()
pd.pm10 <- c()

for (i in 1:nrow(data)) {
  
  print(i)
  tmp <- data[which(data$DOY == (data$DOY[i]-1)),]
  pd.pm <- c(pd.pm, tmp$PM[1])
  pd.pm10 <- c(pd.pm10, tmp$PM10[1])
  
}

data$PM_Lag <- pd.pm
data$PM10_Lag <- pd.pm10

# Main analysis with full sample

data$TICKETED2 <- data$TICKETED + data$ARRESTED # All citations

# Save this version of data to minimize the intense crying

write.csv(data, paste0(direc, 'data/final_data.csv'), row.names = FALSE)

# OLS REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket.xxx <- lm(TICKETED2 ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

ticket.x10 <- lm(TICKETED2 ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xticket.xxx <- coeftest(ticket.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket.x10 <- coeftest(ticket.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested.xxx <- lm(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

arrested.x10 <- lm(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xarrested.xxx <- coeftest(arrested.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested.x10 <- coeftest(arrested.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx.xxx <- lm(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

arrestedx.x10 <- lm(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

xarrestedx.xxx <- coeftest(arrestedx.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx.x10 <- coeftest(arrestedx.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched.xxx <- lm(SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

searched.x10 <- lm(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xsearched.xxx <- coeftest(searched.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched.x10 <- coeftest(searched.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle.xxx <- lm(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

vehicle.x10 <- lm(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

xvehicle.xxx <- coeftest(vehicle.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle.x10 <- coeftest(vehicle.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver.xxx <- lm(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

driver.x10 <- lm(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xdriver.xxx <- coeftest(driver.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver.x10 <- coeftest(driver.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force.xxx <- lm(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

force.x10 <- lm(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

xforce.xxx <- coeftest(force.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce.x10 <- coeftest(force.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex.xxx <- lm(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

forcex.x10 <- lm(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

xforcex.xxx <- coeftest(forcex.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex.x10 <- coeftest(forcex.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# IV REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED2 ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

ticket_iv.x10 <- ivreg(TICKETED2 ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

searched_iv.x10 <- ivreg(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

force_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
          xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv.txt'), row.names = FALSE)

write.csv(stargazer(xticket.xxx, xarrested.xxx, xarrestedx.xxx, xsearched.xxx, xvehicle.xxx, xdriver.xxx, xforce.xxx, xforcex.xxx,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
          force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj.txt'), row.names = FALSE)

write.csv(stargazer(ticket.xxx, arrested.xxx, arrestedx.xxx, searched.xxx, vehicle.xxx, driver.xxx, force.xxx, forcex.xxx,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
          xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv.txt'), row.names = FALSE)

write.csv(stargazer(xticket.x10, xarrested.x10, xarrestedx.x10, xsearched.x10, xvehicle.x10, xdriver.x10, xforce.x10, xforcex.x10,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
          force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj.txt'), row.names = FALSE)

write.csv(stargazer(ticket.x10, arrested.x10, arrestedx.x10, searched.x10, vehicle.x10, driver.x10, force.x10, forcex.x10,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_unadj.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                 sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                               c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                 sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats) <- c('PM', 'PM10')

write.csv(f_stats, paste0(direc, 'results/f_stats.txt'), row.names = FALSE)

# IV robustness test omitting XXX as a second instrument

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED2 ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY), data = data)

ticket_iv.x10 <- ivreg(TICKETED2 ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY), data = data)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY), data = data)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY), data = data)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY), data = data[which(data$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY), data = data[which(data$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY), data = data)

searched_iv.x10 <- ivreg(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY), data = data)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY), data = data)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY), data = data)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY), data = data)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY), data = data)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY), data = data)

force_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY), data = data)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY), data = data[which(data$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY), data = data[which(data$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_robust_no_xxx.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_robust_no_xxx.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_robust_no_xxx.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_robust_no_xxx.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                 sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                               c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                 sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats) <- c('PM', 'PM10')

write.csv(f_stats, paste0(direc, 'results/f_stats_robust_no_xxx.txt'), row.names = FALSE)

# Another IV robustness test (restricted sample in time)

data2 <- data[which(data$DAY %in% c(136:205)),]

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED2 ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data2)

ticket_iv.x10 <- ivreg(TICKETED2 ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data2)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data2)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data2)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data2[which(data2$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data2[which(data2$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data2)

searched_iv.x10 <- ivreg(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data2)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data2)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data2)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data2)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data2)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data2)

force_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data2)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data2[which(data2$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data2[which(data2$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_robust_restricted.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_robust_restricted.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_robust_restricted.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_robust_restricted.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                 sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7]),
                               c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                 sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7])))

colnames(f_stats) <- c('PM', 'PM10')

write.csv(f_stats, paste0(direc, 'results/f_stats_robust_restricted.txt'), row.names = FALSE)

# Linear PM robustness test

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED2 ~ PM + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - PM - I(PM^2) - I(PM^3) + PM_PY + I(PM_PY^2) + XXX, data = data)

ticket_iv.x10 <- ivreg(TICKETED2 ~ PM10 + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - PM10 - I(PM10^2) - I(PM10^3) + PM10_PY + I(PM10_PY^2) + XXX, data = data)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ PM + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - PM - I(PM^2) - I(PM^3) + PM_PY + I(PM_PY^2) + XXX, data = data)

arrested_iv.x10 <- ivreg(ARRESTED ~ PM10 + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - PM10 - I(PM10^2) - I(PM10^3) + PM10_PY + I(PM10_PY^2) + XXX, data = data)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ PM + log(PM_Lag) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - PM - I(PM^2) - I(PM^3) + PM_PY + I(PM_PY^2) + XXX, data = data[which(data$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ PM10 + log(PM10_Lag) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - PM10 - I(PM10^2) - I(PM10^3) + PM10_PY + I(PM10_PY^2) + XXX, data = data[which(data$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ PM + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - PM - I(PM^2) - I(PM^3) + PM_PY + I(PM_PY^2) + XXX, data = data)

searched_iv.x10 <- ivreg(SEARCHED ~ PM10 + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - PM10 - I(PM10^2) - I(PM10^3) + PM10_PY + I(PM10_PY^2) + XXX, data = data)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ PM + log(PM_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - PM - I(PM^2) - I(PM^3) + PM_PY + I(PM_PY^2) + XXX, data = data)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ PM10 + log(PM10_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - PM10 - I(PM10^2) - I(PM10^3) + PM10_PY + I(PM10_PY^2) + XXX, data = data)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ PM + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - PM - I(PM^2) - I(PM^3) + PM_PY + I(PM_PY^2) + XXX, data = data)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ PM10 + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - PM10 - I(PM10^2) - I(PM10^3) + PM10_PY + I(PM10_PY^2) + XXX, data = data)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ PM + log(PM_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - PM - I(PM^2) - I(PM^3) + PM_PY + I(PM_PY^2) + XXX, data = data)

force_iv.x10 <- ivreg(FORCE ~ PM10 + log(PM10_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - PM10 - I(PM10^2) - I(PM10^3) + PM10_PY + I(PM10_PY^2) + XXX, data = data)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ PM + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - PM - I(PM^2) - I(PM^3) + PM_PY + I(PM_PY^2) + XXX, data = data[which(data$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ PM10 + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - PM10 - I(PM10^2) - I(PM10^3) + PM10_PY + I(PM10_PY^2) + XXX, data = data[which(data$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_robust_linear.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_robust_linear.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_robust_linear.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_robust_linear.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                 sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                               c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                 sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats) <- c('PM', 'PM10')

write.csv(f_stats, paste0(direc, 'results/f_stats_robust_linear.txt'), row.names = FALSE)

# Add weather data to counts

counts <- left_join(counts, data[,which(colnames(data) %in% c('STOP_DATE', 'JURISDICTION', 'Ozone', 'CO', 'TEMP', 'PRCP', 'PM_PY', 'PD', 'WEEK', 'DAY', 'Ozone_PY', 'CO_PY'))], by = c('DATE2' = 'STOP_DATE', 'Jurisdiction' = 'JURISDICTION'))
counts <- counts[!duplicated(counts),]

# Counts analyses

c.lags <- c()
c.lags10 <- c()
c.xxx <- c()

for (i in 1:nrow(counts)) {
  
  print(i)
  
  tmp <- data %>% filter(DOY == counts$DAY[i])
  tmp <- tmp %>% filter(JURISDICTION == counts$Jurisdiction[i])
  
  c.lags <- c(c.lags, tmp$PM_Lag[1])
  c.lags10 <- c(c.lags10, tmp$PM10_Lag[1])
  c.xxx <- c(c.xxx, tmp$XXX[1])
  
}

counts$PM_Lag <- c.lags
counts$PM10_Lag <- c.lags10
counts$XXX <- c.xxx
counts$WEEK <- counts$WEEK.x
counts$DOW <- counts$DAY %% 7

# Stops

ivstops1 <- ivreg(log(Stops+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Agency) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops2 <- ivreg(log(Stops+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Jurisdiction) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops3 <- ivreg(log(Stops+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Agency) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops4 <- ivreg(log(Stops+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Jurisdiction) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops5 <- ivreg(log(Stops+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Agency) + factor(Jurisdiction) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops6 <- ivreg(log(Stops+1) ~ log(PM) + log(Ozone) + log(CO) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Agency) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops7 <- ivreg(log(Stops+1) ~ log(PM) + log(Ozone) + log(CO) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Jurisdiction) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops8 <- ivreg(log(Stops+1) ~ log(PM) + log(Ozone) + log(CO) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Agency) | . - log(PM) + log(PM_PY) + XXX, data = counts)

ivstops1x <- coeftest(ivstops1, vcov = vcovCL(ivstops1, type = 'HC1'))
ivstops2x <- coeftest(ivstops2, vcov = vcovCL(ivstops2, type = 'HC1'))
ivstops3x <- coeftest(ivstops3, vcov = vcovCL(ivstops3, type = 'HC1'))
ivstops4x <- coeftest(ivstops4, vcov = vcovCL(ivstops4, type = 'HC1'))
ivstops5x <- coeftest(ivstops5, vcov = vcovCL(ivstops5, type = 'HC1'))
ivstops6x <- coeftest(ivstops6, vcov = vcovCL(ivstops6, type = 'HC1'))
ivstops7x <- coeftest(ivstops7, vcov = vcovCL(ivstops7, type = 'HC1'))
ivstops8x <- coeftest(ivstops8, vcov = vcovCL(ivstops8, type = 'HC1'))

write.csv(stargazer(ivstops1, ivstops2, ivstops3, ivstops4, ivstops5, ivstops6, ivstops7, ivstops8, type = 'text',
                    omit = c('WEEK', 'DOW', 'Agency', 'Jurisdiction')), paste0(direc, 'results/pm_iv_counts_unadj.txt'), row.names = FALSE)

write.csv(stargazer(ivstops1x, ivstops2x, ivstops3x, ivstops4x, ivstops5x, ivstops6x, ivstops7x, ivstops8x, type = 'text',
                    omit = c('WEEK', 'DOW', 'Agency', 'Jurisdiction')), paste0(direc, 'results/pm_iv_counts.txt'), row.names = FALSE)

sum_stops1 <- summary(ivstops1, diagnostics = TRUE)
sum_stops2 <- summary(ivstops2, diagnostics = TRUE)
sum_stops3 <- summary(ivstops3, diagnostics = TRUE)
sum_stops4 <- summary(ivstops4, diagnostics = TRUE)
sum_stops5 <- summary(ivstops5, diagnostics = TRUE)

sum_stops6 <- summary(ivstops6, diagnostics = TRUE)
sum_stops7 <- summary(ivstops7, diagnostics = TRUE)
sum_stops8 <- summary(ivstops8, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_stops1[[12]][7], sum_stops2[[12]][7], sum_stops3[[12]][7], sum_stops4[[12]][7],
                                 sum_stops5[[12]][7], sum_stops6[[12]][7], sum_stops7[[12]][7], sum_stops8[[12]][7])))

colnames(f_stats) <- c('PM')

write.csv(f_stats, paste0(direc, 'results/f_stats_stops.txt'), row.names = FALSE)

# Re-run the IV analysis for same-county residents only

samesies <- data %>% filter(RESIDENCY == 'RESIDENT OF CITY/COUNTY OF STOP')

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesies)

ticket_iv.x10 <- ivreg(TICKETED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesies)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesies)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesies)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + log(Searched+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesies[which(samesies$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + log(Searched+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesies[which(samesies$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesies)

searched_iv.x10 <- ivreg(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesies)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesies)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesies)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesies)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesies)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesies)

force_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesies)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + log(Arrests+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesies[which(samesies$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + log(Arrests+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesies[which(samesies$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_local.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_local.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_local.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_local.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats3 <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                  sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                                c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                  sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats3) <- c('PM', 'PM10')

write.csv(f_stats3, paste0(direc, 'results/f_stats_local.txt'), row.names = FALSE)

# Re-run the IV analysis for non-same-county residents only

samesiesish <- data %>% filter(RESIDENCY != 'RESIDENT OF CITY/COUNTY OF STOP')

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesish)

ticket_iv.x10 <- ivreg(TICKETED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesish)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesish)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesish)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + log(Searched+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesish[which(samesiesish$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + log(Searched+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesish[which(samesiesish$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesish)

searched_iv.x10 <- ivreg(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesish)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesish)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesish)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesish)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesish)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesish)

force_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesish)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + log(Arrests+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesish[which(samesiesish$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + log(Arrests+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesish[which(samesiesish$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_non_local.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_non_local.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_non_local.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_non_local.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats3 <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                  sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                                c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                  sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats3) <- c('PM', 'PM10')

write.csv(f_stats3, paste0(direc, 'results/f_stats_non_local.txt'), row.names = FALSE)

# Re-run the IV analysis for out-of-state drivers only

samesiesishbutnot <- data %>% filter(RESIDENCY == 'OUT OF STATE RESIDENT')

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesishbutnot)

ticket_iv.x10 <- ivreg(TICKETED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesishbutnot)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesishbutnot)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesishbutnot)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Stops+1) + log(Searched+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesishbutnot[which(samesiesishbutnot$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + log(Searched+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesishbutnot[which(samesiesishbutnot$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesishbutnot)

searched_iv.x10 <- ivreg(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesishbutnot)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesishbutnot)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesishbutnot)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesishbutnot)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesishbutnot)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesishbutnot)

force_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesishbutnot)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Stops+1) + log(Arrests+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = samesiesishbutnot[which(samesiesishbutnot$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Stops+1) + log(Arrests+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = samesiesishbutnot[which(samesiesishbutnot$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_out_of_state.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_out_of_state.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_out_of_state.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_out_of_state.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats3 <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                  sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                                c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                  sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats3) <- c('PM', 'PM10')

write.csv(f_stats3, paste0(direc, 'results/f_stats_out_of_state.txt'), row.names = FALSE)

# Re-run main analysis (OLS and IV) with all pollutants (includes Ozone and CO in addition to PM)

ticket.xxx <- lm(TICKETED2 ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

ticket.x10 <- lm(TICKETED2 ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xticket.xxx <- coeftest(ticket.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket.x10 <- coeftest(ticket.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested.xxx <- lm(ARRESTED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

arrested.x10 <- lm(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xarrested.xxx <- coeftest(arrested.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested.x10 <- coeftest(arrested.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx.xxx <- lm(ARRESTED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

arrestedx.x10 <- lm(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

xarrestedx.xxx <- coeftest(arrestedx.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx.x10 <- coeftest(arrestedx.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched.xxx <- lm(SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

searched.x10 <- lm(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xsearched.xxx <- coeftest(searched.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched.x10 <- coeftest(searched.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle.xxx <- lm(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

vehicle.x10 <- lm(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

xvehicle.xxx <- coeftest(vehicle.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle.x10 <- coeftest(vehicle.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver.xxx <- lm(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

driver.x10 <- lm(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xdriver.xxx <- coeftest(driver.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver.x10 <- coeftest(driver.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force.xxx <- lm(FORCE ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

force.x10 <- lm(FORCE ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

xforce.xxx <- coeftest(force.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce.x10 <- coeftest(force.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex.xxx <- lm(FORCE ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

forcex.x10 <- lm(FORCE ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

xforcex.xxx <- coeftest(forcex.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex.x10 <- coeftest(forcex.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# IV REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED2 ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

ticket_iv.x10 <- ivreg(TICKETED2 ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

searched_iv.x10 <- ivreg(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

force_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_all_pollutants.txt'), row.names = FALSE)

write.csv(stargazer(xticket.xxx, xarrested.xxx, xarrestedx.xxx, xsearched.xxx, xvehicle.xxx, xdriver.xxx, xforce.xxx, xforcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_all_pollutants.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_all_pollutants.txt'), row.names = FALSE)

write.csv(stargazer(ticket.xxx, arrested.xxx, arrestedx.xxx, searched.xxx, vehicle.xxx, driver.xxx, force.xxx, forcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_all_pollutants.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_all_pollutants.txt'), row.names = FALSE)

write.csv(stargazer(xticket.x10, xarrested.x10, xarrestedx.x10, xsearched.x10, xvehicle.x10, xdriver.x10, xforce.x10, xforcex.x10,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_all_pollutants.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_all_pollutants.txt'), row.names = FALSE)

write.csv(stargazer(ticket.x10, arrested.x10, arrestedx.x10, searched.x10, vehicle.x10, driver.x10, force.x10, forcex.x10,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_unadj_all_pollutants.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                 sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                               c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                 sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats) <- c('PM', 'PM10')

write.csv(f_stats, paste0(direc, 'results/f_stats_all_pollutants.txt'), row.names = FALSE)

# EVERYTHING MODELS !!!

# OLS REGRESSIONS

ticket.xxx <- lm(TICKETED2 ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

ticket.x10 <- lm(TICKETED2 ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xticket.xxx <- coeftest(ticket.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket.x10 <- coeftest(ticket.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested.xxx <- lm(ARRESTED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

arrested.x10 <- lm(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xarrested.xxx <- coeftest(arrested.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested.x10 <- coeftest(arrested.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx.xxx <- lm(ARRESTED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

arrestedx.x10 <- lm(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

xarrestedx.xxx <- coeftest(arrestedx.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx.x10 <- coeftest(arrestedx.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched.xxx <- lm(SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

searched.x10 <- lm(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xsearched.xxx <- coeftest(searched.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched.x10 <- coeftest(searched.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle.xxx <- lm(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

vehicle.x10 <- lm(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

xvehicle.xxx <- coeftest(vehicle.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle.x10 <- coeftest(vehicle.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver.xxx <- lm(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

driver.x10 <- lm(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xdriver.xxx <- coeftest(driver.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver.x10 <- coeftest(driver.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force.xxx <- lm(FORCE ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

force.x10 <- lm(FORCE ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

xforce.xxx <- coeftest(force.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce.x10 <- coeftest(force.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex.xxx <- lm(FORCE ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

forcex.x10 <- lm(FORCE ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

xforcex.xxx <- coeftest(forcex.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex.x10 <- coeftest(forcex.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# IV REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED2 ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

ticket_iv.x10 <- ivreg(TICKETED2 ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

searched_iv.x10 <- ivreg(SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

force_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_everything.txt'), row.names = FALSE)

write.csv(stargazer(xticket.xxx, xarrested.xxx, xarrestedx.xxx, xsearched.xxx, xvehicle.xxx, xdriver.xxx, xforce.xxx, xforcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_everything.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_everything.txt'), row.names = FALSE)

write.csv(stargazer(ticket.xxx, arrested.xxx, arrestedx.xxx, searched.xxx, vehicle.xxx, driver.xxx, force.xxx, forcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_everything.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_everything.txt'), row.names = FALSE)

write.csv(stargazer(xticket.x10, xarrested.x10, xarrestedx.x10, xsearched.x10, xvehicle.x10, xdriver.x10, xforce.x10, xforcex.x10,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_everything.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_everything.txt'), row.names = FALSE)

write.csv(stargazer(ticket.x10, arrested.x10, arrestedx.x10, searched.x10, vehicle.x10, driver.x10, force.x10, forcex.x10,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_unadj_everything.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                 sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                               c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                 sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats) <- c('PM', 'PM10')

write.csv(f_stats, paste0(direc, 'results/f_stats_everything.txt'), row.names = FALSE)

# Figure for stops by day to compare to pollution levels by day

day <- c()
ds <- c()

for (i in 1:max(data$DAY)) {
  
  tmp <- data %>% filter(DAY == i)
  day <- c(day, i)
  ds <- c(ds, nrow(tmp))
  
}

fdfs <- as.data.frame(cbind(day, ds))
colnames(fdfs) <- c('Day', 'Stops')

mark <- lm(Stops ~ 1, data = fdfs)

msms <- msmFit(mark, k = 2, sw = rep(TRUE, 2))

plotProb(msms, which = 1)

fdfs$Color <- c(rep('Baseline Pollution', 135), rep('Canadian Wildfires', 70), rep('Baseline Pollution', 160))

ggplot(data = fdfs, aes(x = Day, y = Stops)) + 
  theme_bw() + 
  ggtitle('Daily Traffic Stop Counts') + 
  ylab('Count') +
  xlab('Day') +
  geom_path(aes(x = Day, y = Stops, group = 1, color = Color), size = 1, alpha = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = 136) +
  geom_vline(xintercept = 206) +
  ylim(c(0,7000)) +
  scale_x_continuous(breaks = seq(0, 365, 30), labels = seq(0, 365, 30)) +
  scale_color_manual(breaks = c('Baseline Pollution', 'Canadian Wildfires'), values = c('orange2', 'red4'), name = 'Regime') +
  annotate('text', x = 171, y = 6200, label = 'Canadian Wildfires') +
  annotate('text', x = 67, y = 6200, label = 'Baseline Pollution') +
  annotate('text', x = 287, y = 6200, label = 'Baseline Pollution')

# Make maps of counties included in sample for PM and by annual stops with opacity zero

va.co <- counties(state = 'VA')

yep <- c('Fairfax County', 'Loudoun County', 'Frederick County', 'Virginia Beach city',
         'Hampton city', 'Norfolk city', 'Rockingham County', 'Salem city',
         'Henrico County', 'Richmond city', 'Chesterfield County',
         'Charles City County', 'Lynchburg city', 'Roanoke County')

va.co$Color <- ifelse(va.co$NAMELSAD %in% yep, 1, 0)

leaflet(va.co$geometry) %>% addTiles() %>% addPolygons(weight = 0.4, smoothFactor = 1.0, opacity = 1.0, fillOpacity = va.co$Color, color = 'black', fillColor = 'red')

# Creating summary statistics tables

#data$Female <- ifelse(data$GENDER == 'FEMALE', 1, 0)
#data$Male <- ifelse(data$GENDER == 'FEMALE', 0, 1)
#data$Hispanic <- ifelse(data$ETHNICITY == 'HISPANIC OR LATINO', 1, 0)
#data$English <- ifelse(data$ENGLISH.SPEAKING == 'YES', 1, 0)
#data$Asian <- ifelse(data$RACE == 'ASIAN OR NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 1, 0)
#data$Black <- ifelse(data$RACE == 'BLACK OR AFRICAN AMERICAN', 1, 0)
#data$Native <- ifelse(data$RACE == 'AMERICAN INDIAN OR ALASKA NATIVE', 1, 0)
#data$White <- ifelse(data$RACE == 'WHITE', 1, 0)

#dfsum <- data[,c(104, 68, 66, 64, 65, 70, 52, 53, 57, 58, 54, 56, 27)]

#colnames(dfsum) <- c('Ticketed', 'Arrested', 'Searched', 'Vehicle Searched',
#                     'Person Searched', 'Force Used', 'Particulate Matter (PM2.5)',
#                     'Particulate Matter (PM10)', 'PY Particulate Matter (PM2.5)',
#                     'PY Particulate Matter (PM10)', 'Ozone', 'Carbon Monoxide', 'Temperature')

#dfsum2 <- dfsum[which(complete.cases(dfsum)),]

#datasummary_skim(dfsum2, fmt = '%.3f')

# Show that in sample counties are no different from all of VA and rest of VA in a table

dfsum$White <- data$White
dfsum$Black <- data$Black
dfsum$Asian <- data$Asian
dfsum$Hispanic <- data$Hispanic
dfsum$Female <- data$Female

dfsum2 <- dfsum[which(complete.cases(dfsum)),]
dfsum3 <- dfsum[-which(complete.cases(dfsum)),]

v0 <- c('Ticketed', '', 'Arrested', '', 'Searched', '', 'Person Searched', '',
        'Vehicle Searched', '', 'Force Used', '', 'Temperature', '',
        'White', '', 'Black', '', 'Asian', '', 'Hispanic', '', 'Female', '')

v1 <- c(mean(dfsum$Ticketed, na.rm = TRUE), sd(dfsum$Ticketed, na.rm = TRUE), 
        mean(dfsum$Arrested, na.rm = TRUE), sd(dfsum$Arrested, na.rm = TRUE), 
        mean(dfsum$Searched, na.rm = TRUE), sd(dfsum$Searched, na.rm = TRUE), 
        mean(dfsum$`Person Searched`, na.rm = TRUE), sd(dfsum$`Person Searched`, na.rm = TRUE), 
        mean(dfsum$`Vehicle Searched`, na.rm = TRUE), sd(dfsum$`Vehicle Searched`, na.rm = TRUE), 
        mean(dfsum$`Force Used`, na.rm = TRUE), sd(dfsum$`Force Used`, na.rm = TRUE), 
        mean(dfsum$Temperature, na.rm = TRUE), sd(dfsum$Temperature, na.rm = TRUE), 
        mean(dfsum$White, na.rm = TRUE), sd(dfsum$White, na.rm = TRUE), 
        mean(dfsum$Black, na.rm = TRUE), sd(dfsum$Black, na.rm = TRUE), 
        mean(dfsum$Asian, na.rm = TRUE), sd(dfsum$Asian, na.rm = TRUE), 
        mean(dfsum$Hispanic, na.rm = TRUE), sd(dfsum$Hispanic, na.rm = TRUE), 
        mean(dfsum$Female, na.rm = TRUE), sd(dfsum$Female, na.rm = TRUE))

v2 <- c(mean(dfsum2$Ticketed, na.rm = TRUE), sd(dfsum2$Ticketed, na.rm = TRUE), 
        mean(dfsum2$Arrested, na.rm = TRUE), sd(dfsum2$Arrested, na.rm = TRUE), 
        mean(dfsum2$Searched, na.rm = TRUE), sd(dfsum2$Searched, na.rm = TRUE), 
        mean(dfsum2$`Person Searched`, na.rm = TRUE), sd(dfsum2$`Person Searched`, na.rm = TRUE), 
        mean(dfsum2$`Vehicle Searched`, na.rm = TRUE), sd(dfsum2$`Vehicle Searched`, na.rm = TRUE), 
        mean(dfsum2$`Force Used`, na.rm = TRUE), sd(dfsum2$`Force Used`, na.rm = TRUE), 
        mean(dfsum2$Temperature, na.rm = TRUE), sd(dfsum2$Temperature, na.rm = TRUE), 
        mean(dfsum2$White, na.rm = TRUE), sd(dfsum2$White, na.rm = TRUE), 
        mean(dfsum2$Black, na.rm = TRUE), sd(dfsum2$Black, na.rm = TRUE), 
        mean(dfsum2$Asian, na.rm = TRUE), sd(dfsum2$Asian, na.rm = TRUE), 
        mean(dfsum2$Hispanic, na.rm = TRUE), sd(dfsum2$Hispanic, na.rm = TRUE), 
        mean(dfsum2$Female, na.rm = TRUE), sd(dfsum2$Female, na.rm = TRUE))

v3 <- c(mean(dfsum3$Ticketed, na.rm = TRUE), sd(dfsum3$Ticketed, na.rm = TRUE), 
        mean(dfsum3$Arrested, na.rm = TRUE), sd(dfsum3$Arrested, na.rm = TRUE), 
        mean(dfsum3$Searched, na.rm = TRUE), sd(dfsum3$Searched, na.rm = TRUE), 
        mean(dfsum3$`Person Searched`, na.rm = TRUE), sd(dfsum3$`Person Searched`, na.rm = TRUE), 
        mean(dfsum3$`Vehicle Searched`, na.rm = TRUE), sd(dfsum3$`Vehicle Searched`, na.rm = TRUE), 
        mean(dfsum3$`Force Used`, na.rm = TRUE), sd(dfsum3$`Force Used`, na.rm = TRUE), 
        mean(dfsum3$Temperature, na.rm = TRUE), sd(dfsum3$Temperature, na.rm = TRUE), 
        mean(dfsum3$White, na.rm = TRUE), sd(dfsum3$White, na.rm = TRUE), 
        mean(dfsum3$Black, na.rm = TRUE), sd(dfsum3$Black, na.rm = TRUE), 
        mean(dfsum3$Asian, na.rm = TRUE), sd(dfsum3$Asian, na.rm = TRUE), 
        mean(dfsum3$Hispanic, na.rm = TRUE), sd(dfsum3$Hispanic, na.rm = TRUE), 
        mean(dfsum3$Female, na.rm = TRUE), sd(dfsum3$Female, na.rm = TRUE))

compdf <- as.data.frame(cbind(v1, v2, v3))
colnames(compdf) <- c('Virginia', 'Sample', 'Difference')

compdf$Virginia <- round(compdf$Virginia, 3)
compdf$Sample <- round(compdf$Sample, 3)
compdf$Difference <- round(compdf$Difference, 3)

compdf <- cbind(v0, compdf)
colnames(compdf)[1] <- 'Variable'

write.csv(compdf, paste0(direc, 'results/comparison_statistics.txt'), row.names = FALSE)

# OK, finally time for the crashes

cds <- c()

for (i in 1:nrow(crashes)) {
  
  print(i)
  cds <- c(cds, strsplit(crashes$Crash.Date[i], ' ')[[1]][1])
  
}

crashes$CD <- cd
cdid <- c()

for (i in 1:nrow(crashes)) {
  
  print(i)
  poop <- strsplit(crashes$CD[i], '/')
  cdid <- c(cdid, ((100*as.integer(poop[[1]][1])) + as.integer(poop[[1]][2])))
  
}

crashes$CDID <- cdid
unidays <- unique(cdid)[order(unique(cdid))]
cday <- c()

for (i in 1:nrow(crashes)) {
  
  print(i)
  cday <- c(cday, which(unidays == crashes$CDID[i]))
  
}

crashes$Day <- cday

yep2 <- c('029. Fairfax County', '053. Loudoun County', '034. Frederick County', '134. City of Virginia Beach',
          '114. City of Hampton', '122. City of Norfolk', '082. Rockingham County', '129. City of Salem',
          '043. Henrico County', '127. City of Richmond', '020. Chesterfield County',
          '018. Charles City County', '118. City of Lynchburg', '080. Roanoke County')

crashes <- crashes %>% filter(Physical.Juris.Name %in% yep2)

# 51373 / 127596 = 40.26% of crashes in Virginia in 2023 is in this sample

crash.days <- c()

for (name in unique(crashes$Physical.Juris.Name)) {
  
  for (i in 1:365) {
    
    crash.days <- c(crash.days, paste0(name, i))
    
  }
  
}

crashes$JD <- paste0(crashes$Physical.Juris.Name, crashes$Day)

crashes$Young <- ifelse(crashes$Young. == 'Yes', 1, 0)
crashes$Senior <- ifelse(crashes$Senior. == 'Yes', 1, 0)
crashes$Alcohol <- ifelse(crashes$Alcohol. == 'Yes', 1, 0)
crashes$Animal <- ifelse(crashes$Animal.Related. == 'Yes', 1, 0)
crashes$Unrestrained <- ifelse(crashes$Unrestrained. == 'Belted', 0, 1)
crashes$School <- ifelse(crashes$School.Zone == '3. No', 0, 1)
crashes$School2 <- ifelse(crashes$School.Zone == '2. Yes - With School Activity', 1, 0)
crashes$Work <- ifelse(crashes$Work.Zone.Related == '1. Yes', 1, 0)
crashes$Distracted <- ifelse(crashes$Distracted. == '1. Yes', 1, 0)
crashes$Drowsy <- ifelse(crashes$Drowsy. == '1. Yes', 1, 0)
crashes$Drugs <- ifelse(crashes$Drug.Related. == '1. Yes', 1, 0)
crashes$Pedestrian <- ifelse(crashes$Pedestrian. == 'Yes', 1, 0)
crashes$HitRun <- ifelse(crashes$Hitrun. == 'Yes', 1, 0)
crashes$Night <- ifelse(crashes$Night. == 'Yes', 1, 0)

crash.count <- c()
young.count <- c()
senior.count <- c()
alcohol.count <- c()
animal.count <- c()
unrestrained.count <- c()
school.count <- c()
school2.count <- c()
work.count <- c()
distracted.count <- c()
drowsy.count <- c()
drugs.count <- c()
pedestrian.count <- c()
hitrun.count <- c()
night.count <- c()

for (i in 1:length(crash.days)) {
  
  print(i)
  
  tmp <- crashes %>% filter(JD == crash.days[i])
  
  crash.count <- c(crash.count, nrow(tmp))
  young.count <- c(young.count, sum(tmp$Young, na.rm = TRUE))
  senior.count <- c(senior.count, sum(tmp$Senior, na.rm = TRUE))
  alcohol.count <- c(alcohol.count, sum(tmp$Alcohol, na.rm = TRUE))
  animal.count <- c(animal.count, sum(tmp$Animal, na.rm = TRUE))
  unrestrained.count <- c(unrestrained.count, sum(tmp$Unrestrained, na.rm = TRUE))
  school.count <- c(school.count, sum(tmp$School, na.rm = TRUE))
  school2.count <- c(school2.count, sum(tmp$School2, na.rm = TRUE))
  work.count <- c(work.count, sum(tmp$Work, na.rm = TRUE))
  distracted.count <- c(distracted.count, sum(tmp$Distracted, na.rm = TRUE))
  drowsy.count <- c(drowsy.count, sum(tmp$Drowsy, na.rm = TRUE))
  drugs.count <- c(drugs.count, sum(tmp$Drugs, na.rm = TRUE))
  pedestrian.count <- c(pedestrian.count, sum(tmp$Pedestrian, na.rm = TRUE))
  hitrun.count <- c(hitrun.count, sum(tmp$HitRun, na.rm = TRUE))
  night.count <- c(night.count, sum(tmp$Night, na.rm = TRUE))
  
}

cdf <- as.data.frame(cbind(crash.count, young.count, senior.count, alcohol.count, animal.count,
                           unrestrained.count, school.count, school2.count, work.count, distracted.count,
                           drowsy.count, drugs.count, pedestrian.count, hitrun.count, night.count))

yep3 <- c('FAIRFAX CO', 'LOUDOUN CO', 'FREDERICK CO', 'VIRGINIA BEACH',
          'HAMPTON', 'NORFOLK', 'ROCKINGHAM CO', 'SALEM',
          'HENRICO CO', 'RICHMOND CITY', 'CHESTERFIELD CO',
          'CHARLES CITY CO', 'LYNCHBURG', 'ROANOKE CO')

old.name <- c()

for (i in 1:nrow(crashes)) {
  
  print(i)
  old.name <- c(old.name, yep3[which(yep2 == crashes$Physical.Juris.Name[i])])
  
}

crashes$OLD <- old.name

crash.days2 <- c()

for (name in unique(crashes$OLD)) {
  
  for (i in 1:365) {
    
    crash.days2 <- c(crash.days2, paste0(name, i))
    
  }
  
}

data$CD2 <- paste0(data$JURISDICTION, data$DAY)

cpm <- c()
cw <- c()
cdow <- c()
cjur <- c()
cpy <- c()
cxxx <- c()
clag <- c()

for (i in 1:length(crash.days2)) {
  
  print(i)
  tmp <- data %>% filter(CD2 == crash.days2[i])
  cpm <- c(cpm, tmp$PM[1])
  cw <- c(cw, tmp$WEEK[1])
  cdow <- c(cdow, tmp$DOW[1])
  cjur <- c(cjur, tmp$JURISDICTION[1])
  cpy <- c(cpy, tmp$PM_PY[1])
  cxxx <- c(cxxx, tmp$XXX[1])
  clag <- c(clag, tmp$PM_Lag[1])
  
}

cdf$PM <- cpm
cdf$WEEK <- cw
cdf$DOW <- cdow
cdf$JURISDIC <- cjur
cdf$PM_PY <- cpy
cdf$XXX <- cxxx
cdf$PM_Lag <- clag

# Crashes regressions with lags

c1 <- lm(crash.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c2 <- lm(school.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c3 <- lm(school2.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c4 <- lm(work.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c5 <- lm(unrestrained.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c6 <- lm(alcohol.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c7 <- lm(animal.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c8 <- lm(pedestrian.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c9 <- lm(hitrun.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c10 <- lm(night.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c11 <- lm(young.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c12 <- lm(senior.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)

ic1 <- ivreg(crash.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic2 <- ivreg(school.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic3 <- ivreg(school2.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic4 <- ivreg(work.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic5 <- ivreg(unrestrained.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic6 <- ivreg(alcohol.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic7 <- ivreg(animal.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic8 <- ivreg(pedestrian.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic9 <- ivreg(hitrun.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic10 <- ivreg(night.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic11 <- ivreg(young.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic12 <- ivreg(senior.count ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)

c1x <- coeftest(c1, vcov = vcovCL(c1), type = 'HC1')
c2x <- coeftest(c2, vcov = vcovCL(c2), type = 'HC1')
c3x <- coeftest(c3, vcov = vcovCL(c3), type = 'HC1')
c4x <- coeftest(c4, vcov = vcovCL(c4), type = 'HC1')
c5x <- coeftest(c5, vcov = vcovCL(c5), type = 'HC1')
c6x <- coeftest(c6, vcov = vcovCL(c6), type = 'HC1')
c7x <- coeftest(c7, vcov = vcovCL(c7), type = 'HC1')
c8x <- coeftest(c8, vcov = vcovCL(c8), type = 'HC1')
c9x <- coeftest(c9, vcov = vcovCL(c9), type = 'HC1')
c10x <- coeftest(c10, vcov = vcovCL(c10), type = 'HC1')
c11x <- coeftest(c11, vcov = vcovCL(c11), type = 'HC1')
c12x <- coeftest(c12, vcov = vcovCL(c12), type = 'HC1')

ic1x <- coeftest(ic1, vcov = vcovCL(ic1), type = 'HC1')
ic2x <- coeftest(ic2, vcov = vcovCL(ic2), type = 'HC1')
ic3x <- coeftest(ic3, vcov = vcovCL(ic3), type = 'HC1')
ic4x <- coeftest(ic4, vcov = vcovCL(ic4), type = 'HC1')
ic5x <- coeftest(ic5, vcov = vcovCL(ic5), type = 'HC1')
ic6x <- coeftest(ic6, vcov = vcovCL(ic6), type = 'HC1')
ic7x <- coeftest(ic7, vcov = vcovCL(ic7), type = 'HC1')
ic8x <- coeftest(ic8, vcov = vcovCL(ic8), type = 'HC1')
ic9x <- coeftest(ic9, vcov = vcovCL(ic9), type = 'HC1')
ic10x <- coeftest(ic10, vcov = vcovCL(ic10), type = 'HC1')
ic11x <- coeftest(ic11, vcov = vcovCL(ic11), type = 'HC1')
ic12x <- coeftest(ic12, vcov = vcovCL(ic12), type = 'HC1')

# Saving results

write.csv(stargazer(ic1x, ic2x, ic3x, ic4x, ic5x, ic6x, ic7x, ic8x, ic9x, ic10x, ic11x, ic12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_crashes_lag.txt'), row.names = FALSE)

write.csv(stargazer(c1x, c2x, c3x, c4x, c5x, c6x, c7x, c8x, c9x, c10x, c11x, c12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_crashes_lag.txt'), row.names = FALSE)

write.csv(stargazer(ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, ic9, ic10, ic11, ic12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_unadj_crashes_lag.txt'), row.names = FALSE)

write.csv(stargazer(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_crashes_lag.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum1 <- summary(ic1, diagnostics = TRUE)
sum2 <- summary(ic2, diagnostics = TRUE)
sum3 <- summary(ic3, diagnostics = TRUE)
sum4 <- summary(ic4, diagnostics = TRUE)
sum5 <- summary(ic5, diagnostics = TRUE)
sum6 <- summary(ic6, diagnostics = TRUE)
sum7 <- summary(ic7, diagnostics = TRUE)
sum8 <- summary(ic8, diagnostics = TRUE)
sum9 <- summary(ic9, diagnostics = TRUE)
sum10 <- summary(ic10, diagnostics = TRUE)
sum11 <- summary(ic11, diagnostics = TRUE)
sum12 <- summary(ic12, diagnostics = TRUE)

f_stats <- as.data.frame(c(sum1[[12]][7], sum2[[12]][7], sum3[[12]][7], sum4[[12]][7], sum5[[12]][7], sum6[[12]][7],
                           sum7[[12]][7], sum8[[12]][7], sum9[[12]][7], sum10[[12]][7], sum11[[12]][7], sum12[[12]][7]))

colnames(f_stats) <- c('PM')

write.csv(f_stats, paste0(direc, 'results/f_stats_crashes_lag.txt'), row.names = FALSE)

# Crashes regressions without lags

c1 <- lm(crash.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c2 <- lm(school.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c3 <- lm(school2.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c4 <- lm(work.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c5 <- lm(unrestrained.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c6 <- lm(alcohol.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c7 <- lm(animal.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c8 <- lm(pedestrian.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c9 <- lm(hitrun.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c10 <- lm(night.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c11 <- lm(young.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c12 <- lm(senior.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)

ic1 <- ivreg(crash.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic2 <- ivreg(school.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic3 <- ivreg(school2.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic4 <- ivreg(work.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic5 <- ivreg(unrestrained.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic6 <- ivreg(alcohol.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic7 <- ivreg(animal.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic8 <- ivreg(pedestrian.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic9 <- ivreg(hitrun.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic10 <- ivreg(night.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic11 <- ivreg(young.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic12 <- ivreg(senior.count ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)

c1x <- coeftest(c1, vcov = vcovCL(c1), type = 'HC1')
c2x <- coeftest(c2, vcov = vcovCL(c2), type = 'HC1')
c3x <- coeftest(c3, vcov = vcovCL(c3), type = 'HC1')
c4x <- coeftest(c4, vcov = vcovCL(c4), type = 'HC1')
c5x <- coeftest(c5, vcov = vcovCL(c5), type = 'HC1')
c6x <- coeftest(c6, vcov = vcovCL(c6), type = 'HC1')
c7x <- coeftest(c7, vcov = vcovCL(c7), type = 'HC1')
c8x <- coeftest(c8, vcov = vcovCL(c8), type = 'HC1')
c9x <- coeftest(c9, vcov = vcovCL(c9), type = 'HC1')
c10x <- coeftest(c10, vcov = vcovCL(c10), type = 'HC1')
c11x <- coeftest(c11, vcov = vcovCL(c11), type = 'HC1')
c12x <- coeftest(c12, vcov = vcovCL(c12), type = 'HC1')

ic1x <- coeftest(ic1, vcov = vcovCL(ic1), type = 'HC1')
ic2x <- coeftest(ic2, vcov = vcovCL(ic2), type = 'HC1')
ic3x <- coeftest(ic3, vcov = vcovCL(ic3), type = 'HC1')
ic4x <- coeftest(ic4, vcov = vcovCL(ic4), type = 'HC1')
ic5x <- coeftest(ic5, vcov = vcovCL(ic5), type = 'HC1')
ic6x <- coeftest(ic6, vcov = vcovCL(ic6), type = 'HC1')
ic7x <- coeftest(ic7, vcov = vcovCL(ic7), type = 'HC1')
ic8x <- coeftest(ic8, vcov = vcovCL(ic8), type = 'HC1')
ic9x <- coeftest(ic9, vcov = vcovCL(ic9), type = 'HC1')
ic10x <- coeftest(ic10, vcov = vcovCL(ic10), type = 'HC1')
ic11x <- coeftest(ic11, vcov = vcovCL(ic11), type = 'HC1')
ic12x <- coeftest(ic12, vcov = vcovCL(ic12), type = 'HC1')

# Saving results

write.csv(stargazer(ic1x, ic2x, ic3x, ic4x, ic5x, ic6x, ic7x, ic8x, ic9x, ic10x, ic11x, ic12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_crashes.txt'), row.names = FALSE)

write.csv(stargazer(c1x, c2x, c3x, c4x, c5x, c6x, c7x, c8x, c9x, c10x, c11x, c12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_crashes.txt'), row.names = FALSE)

write.csv(stargazer(ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, ic9, ic10, ic11, ic12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_unadj_crashes.txt'), row.names = FALSE)

write.csv(stargazer(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_crashes.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum1 <- summary(ic1, diagnostics = TRUE)
sum2 <- summary(ic2, diagnostics = TRUE)
sum3 <- summary(ic3, diagnostics = TRUE)
sum4 <- summary(ic4, diagnostics = TRUE)
sum5 <- summary(ic5, diagnostics = TRUE)
sum6 <- summary(ic6, diagnostics = TRUE)
sum7 <- summary(ic7, diagnostics = TRUE)
sum8 <- summary(ic8, diagnostics = TRUE)
sum9 <- summary(ic9, diagnostics = TRUE)
sum10 <- summary(ic10, diagnostics = TRUE)
sum11 <- summary(ic11, diagnostics = TRUE)
sum12 <- summary(ic12, diagnostics = TRUE)

f_stats <- as.data.frame(c(sum1[[12]][7], sum2[[12]][7], sum3[[12]][7], sum4[[12]][7], sum5[[12]][7], sum6[[12]][7],
                           sum7[[12]][7], sum8[[12]][7], sum9[[12]][7], sum10[[12]][7], sum11[[12]][7], sum12[[12]][7]))

colnames(f_stats) <- c('PM')

write.csv(f_stats, paste0(direc, 'results/f_stats_crashes.txt'), row.names = FALSE)

# PLACEBO TESTING

# Placebo tests randomizing outcomes

set.seed(42069)

new.ids <- sample.int(nrow(data), nrow(data))

data$TICKETED2P <- data$TICKETED2[new.ids]
data$ARRESTEDP <- data$ARRESTED[new.ids]
data$SEARCHEDP <- data$SEARCHED[new.ids]
data$VEHICLE_SEARCHEDP <- data$VEHICLE_SEARCHED[new.ids]
data$PERSON_SEARCHEDP <- data$PERSON_SEARCHED[new.ids]
data$FORCEP <- data$FORCE[new.ids]

# OLS REGRESSIONS

ticket.xxx <- lm(TICKETED2P ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

ticket.x10 <- lm(TICKETED2P ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xticket.xxx <- coeftest(ticket.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket.x10 <- coeftest(ticket.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested.xxx <- lm(ARRESTEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

arrested.x10 <- lm(ARRESTEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xarrested.xxx <- coeftest(arrested.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested.x10 <- coeftest(arrested.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx.xxx <- lm(ARRESTEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

arrestedx.x10 <- lm(ARRESTEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

xarrestedx.xxx <- coeftest(arrestedx.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx.x10 <- coeftest(arrestedx.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched.xxx <- lm(SEARCHEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

searched.x10 <- lm(SEARCHEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xsearched.xxx <- coeftest(searched.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched.x10 <- coeftest(searched.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle.xxx <- lm(VEHICLE_SEARCHEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

vehicle.x10 <- lm(VEHICLE_SEARCHEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

xvehicle.xxx <- coeftest(vehicle.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle.x10 <- coeftest(vehicle.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver.xxx <- lm(PERSON_SEARCHEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

driver.x10 <- lm(PERSON_SEARCHEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xdriver.xxx <- coeftest(driver.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver.x10 <- coeftest(driver.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force.xxx <- lm(FORCEP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

force.x10 <- lm(FORCEP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

xforce.xxx <- coeftest(force.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce.x10 <- coeftest(force.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex.xxx <- lm(FORCEP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

forcex.x10 <- lm(FORCEP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

xforcex.xxx <- coeftest(forcex.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex.x10 <- coeftest(forcex.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# IV REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED2P ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

ticket_iv.x10 <- ivreg(TICKETED2P ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

arrested_iv.x10 <- ivreg(ARRESTEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

searched_iv.x10 <- ivreg(SEARCHEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHEDP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

driver_iv.x10 <- ivreg(PERSON_SEARCHEDP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCEP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

force_iv.x10 <- ivreg(FORCEP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCEP ~ log(PM) + log(PM_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCEP ~ log(PM10) + log(PM10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_everything_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(xticket.xxx, xarrested.xxx, xarrestedx.xxx, xsearched.xxx, xvehicle.xxx, xdriver.xxx, xforce.xxx, xforcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_everything_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_everything_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(ticket.xxx, arrested.xxx, arrestedx.xxx, searched.xxx, vehicle.xxx, driver.xxx, force.xxx, forcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_everything_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_everything_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(xticket.x10, xarrested.x10, xarrestedx.x10, xsearched.x10, xvehicle.x10, xdriver.x10, xforce.x10, xforcex.x10,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_everything_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_everything_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(ticket.x10, arrested.x10, arrestedx.x10, searched.x10, vehicle.x10, driver.x10, force.x10, forcex.x10,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_unadj_everything_placebo_outcomes.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                 sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                               c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                 sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats) <- c('PM', 'PM10')

write.csv(f_stats, paste0(direc, 'results/f_stats_everything_placebo_outcomes.txt'), row.names = FALSE)

# Placebo tests randomizing PM

data$Placebo <- data$PM[new.ids]
data$Placebo10 <- data$PM10[new.ids]
data$Placebo_Lag <- data$PM_Lag[new.ids]
data$Placebo10_Lag <- data$PM10_Lag[new.ids]

# OLS REGRESSIONS

ticket.xxx <- lm(TICKETED2 ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

ticket.x10 <- lm(TICKETED2 ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xticket.xxx <- coeftest(ticket.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket.x10 <- coeftest(ticket.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested.xxx <- lm(ARRESTED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

arrested.x10 <- lm(ARRESTED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xarrested.xxx <- coeftest(arrested.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested.x10 <- coeftest(arrested.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx.xxx <- lm(ARRESTED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

arrestedx.x10 <- lm(ARRESTED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                    + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                    + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

xarrestedx.xxx <- coeftest(arrestedx.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx.x10 <- coeftest(arrestedx.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched.xxx <- lm(SEARCHED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

searched.x10 <- lm(SEARCHED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                   + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                   + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xsearched.xxx <- coeftest(searched.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched.x10 <- coeftest(searched.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle.xxx <- lm(VEHICLE_SEARCHED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

vehicle.x10 <- lm(VEHICLE_SEARCHED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                  + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                  + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

xvehicle.xxx <- coeftest(vehicle.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle.x10 <- coeftest(vehicle.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver.xxx <- lm(PERSON_SEARCHED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

driver.x10 <- lm(PERSON_SEARCHED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xdriver.xxx <- coeftest(driver.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver.x10 <- coeftest(driver.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force.xxx <- lm(FORCE ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

force.x10 <- lm(FORCE ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

xforce.xxx <- coeftest(force.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce.x10 <- coeftest(force.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex.xxx <- lm(FORCE ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

forcex.x10 <- lm(FORCE ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                 + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                 + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

xforcex.xxx <- coeftest(forcex.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex.x10 <- coeftest(forcex.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# IV REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED2 ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(Placebo) + log(PM_PY) + XXX, data = data)

ticket_iv.x10 <- ivreg(TICKETED2 ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(Placebo10) + log(PM10_PY) + XXX, data = data)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(Placebo) + log(PM_PY) + XXX, data = data)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(Placebo10) + log(PM10_PY) + XXX, data = data)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(Placebo) + log(PM_PY) + XXX, data = data[which(data$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                          + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                          + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(Placebo10) + log(PM10_PY) + XXX, data = data[which(data$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(Placebo) + log(PM_PY) + XXX, data = data)

searched_iv.x10 <- ivreg(SEARCHED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                         + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                         + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(Placebo10) + log(PM10_PY) + XXX, data = data)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(Placebo) + log(PM_PY) + XXX, data = data)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                        + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                        + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(Placebo10) + log(PM10_PY) + XXX, data = data)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(Placebo) + log(PM_PY) + XXX, data = data)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(Placebo10) + log(PM10_PY) + XXX, data = data)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(Placebo) + log(PM_PY) + XXX, data = data)

force_iv.x10 <- ivreg(FORCE ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                      + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                      + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(Placebo10) + log(PM10_PY) + XXX, data = data)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(Placebo) + log(Placebo_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(Placebo) + log(PM_PY) + XXX, data = data[which(data$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(Placebo10) + log(Placebo10_Lag) + log(Ozone) + log(CO) + log(Stops+1) + MINOR
                       + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + relevel(as.factor(RACE), ref = 'WHITE') + relevel(as.factor(ETHNICITY), ref = 'NOT HISPANIC OR LATINO') + relevel(as.factor(GENDER), ref = 'FEMALE')
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP) + factor(JURISDICTION)
                       + relevel(as.factor(RESIDENCY), ref = 'RESIDENT OF CITY/COUNTY OF STOP') + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(Placebo10) + log(PM10_PY) + XXX, data = data[which(data$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_everything_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(xticket.xxx, xarrested.xxx, xarrestedx.xxx, xsearched.xxx, xvehicle.xxx, xdriver.xxx, xforce.xxx, xforcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_everything_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm_iv_unadj_everything_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(ticket.xxx, arrested.xxx, arrestedx.xxx, searched.xxx, vehicle.xxx, driver.xxx, force.xxx, forcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_everything_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_everything_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(xticket.x10, xarrested.x10, xarrestedx.x10, xsearched.x10, xvehicle.x10, xdriver.x10, xforce.x10, xforcex.x10,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_everything_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'REASON.FOR.STOP')),
          paste0(direc, 'results/pm10_iv_unadj_everything_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(ticket.x10, arrested.x10, arrestedx.x10, searched.x10, vehicle.x10, driver.x10, force.x10, forcex.x10,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_unadj_everything_placebo_pm.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

sum_tix2 <- summary(ticket_iv.x10, diagnostics = TRUE)
sum_arr2 <- summary(arrested_iv.x10, diagnostics = TRUE)
sum_arx2 <- summary(arrestedx_iv.x10, diagnostics = TRUE)
sum_sea2 <- summary(searched_iv.x10, diagnostics = TRUE)
sum_veh2 <- summary(vehicle_iv.x10, diagnostics = TRUE)
sum_dri2 <- summary(driver_iv.x10, diagnostics = TRUE)
sum_for2 <- summary(force_iv.x10, diagnostics = TRUE)
sum_fox2 <- summary(forcex_iv.x10, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                                 sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]),
                               c(sum_tix2[[12]][7], sum_arr2[[12]][7], sum_arx2[[12]][7], sum_sea2[[12]][7],
                                 sum_veh2[[12]][7], sum_dri2[[12]][7], sum_for2[[12]][7], sum_fox2[[12]][7])))

colnames(f_stats) <- c('PM', 'PM10')

write.csv(f_stats, paste0(direc, 'results/f_stats_everything_placebo_pm.txt'), row.names = FALSE)

# Counts analyses

count.new.id <- sample.int(nrow(counts), nrow(counts))
counts$StopsP <- counts$Stops[count.new.id]

# Stops

ivstops1 <- ivreg(log(StopsP+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Agency) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops2 <- ivreg(log(StopsP+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Jurisdiction) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops3 <- ivreg(log(StopsP+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Agency) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops4 <- ivreg(log(StopsP+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Jurisdiction) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops5 <- ivreg(log(StopsP+1) ~ log(PM) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Agency) + factor(Jurisdiction) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops6 <- ivreg(log(StopsP+1) ~ log(PM) + log(Ozone) + log(CO) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Agency) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops7 <- ivreg(log(StopsP+1) ~ log(PM) + log(Ozone) + log(CO) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Jurisdiction) | . - log(PM) + log(PM_PY) + XXX, data = counts)
ivstops8 <- ivreg(log(StopsP+1) ~ log(PM) + log(Ozone) + log(CO) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Agency) | . - log(PM) + log(PM_PY) + XXX, data = counts)

ivstops1x <- coeftest(ivstops1, vcov = vcovCL(ivstops1, type = 'HC1'))
ivstops2x <- coeftest(ivstops2, vcov = vcovCL(ivstops2, type = 'HC1'))
ivstops3x <- coeftest(ivstops3, vcov = vcovCL(ivstops3, type = 'HC1'))
ivstops4x <- coeftest(ivstops4, vcov = vcovCL(ivstops4, type = 'HC1'))
ivstops5x <- coeftest(ivstops5, vcov = vcovCL(ivstops5, type = 'HC1'))
ivstops6x <- coeftest(ivstops6, vcov = vcovCL(ivstops6, type = 'HC1'))
ivstops7x <- coeftest(ivstops7, vcov = vcovCL(ivstops7, type = 'HC1'))
ivstops8x <- coeftest(ivstops8, vcov = vcovCL(ivstops8, type = 'HC1'))

write.csv(stargazer(ivstops1, ivstops2, ivstops3, ivstops4, ivstops5, ivstops6, ivstops7, ivstops8, type = 'text',
                    omit = c('WEEK', 'DOW', 'Agency', 'Jurisdiction')), paste0(direc, 'results/pm_iv_counts_unadj_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(ivstops1x, ivstops2x, ivstops3x, ivstops4x, ivstops5x, ivstops6x, ivstops7x, ivstops8x, type = 'text',
                    omit = c('WEEK', 'DOW', 'Agency', 'Jurisdiction')), paste0(direc, 'results/pm_iv_counts_placebo_outcomes.txt'), row.names = FALSE)

sum_stops1 <- summary(ivstops1, diagnostics = TRUE)
sum_stops2 <- summary(ivstops2, diagnostics = TRUE)
sum_stops3 <- summary(ivstops3, diagnostics = TRUE)
sum_stops4 <- summary(ivstops4, diagnostics = TRUE)
sum_stops5 <- summary(ivstops5, diagnostics = TRUE)

sum_stops6 <- summary(ivstops6, diagnostics = TRUE)
sum_stops7 <- summary(ivstops7, diagnostics = TRUE)
sum_stops8 <- summary(ivstops8, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_stops1[[12]][7], sum_stops2[[12]][7], sum_stops3[[12]][7], sum_stops4[[12]][7],
                                 sum_stops5[[12]][7], sum_stops6[[12]][7], sum_stops7[[12]][7], sum_stops8[[12]][7])))

colnames(f_stats) <- c('PM')

write.csv(f_stats, paste0(direc, 'results/f_stats_stops_placebo_outcomes.txt'), row.names = FALSE)

# Counts analyses

counts$Placebo <- counts$PM[count.new.id]

# Stops

ivstops1 <- ivreg(log(Stops+1) ~ log(Placebo) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Agency) | . - log(Placebo) + log(PM_PY) + XXX, data = counts)
ivstops2 <- ivreg(log(Stops+1) ~ log(Placebo) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Jurisdiction) | . - log(Placebo) + log(PM_PY) + XXX, data = counts)
ivstops3 <- ivreg(log(Stops+1) ~ log(Placebo) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Agency) | . - log(Placebo) + log(PM_PY) + XXX, data = counts)
ivstops4 <- ivreg(log(Stops+1) ~ log(Placebo) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Jurisdiction) | . - log(Placebo) + log(PM_PY) + XXX, data = counts)
ivstops5 <- ivreg(log(Stops+1) ~ log(Placebo) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Agency) + factor(Jurisdiction) | . - log(Placebo) + log(PM_PY) + XXX, data = counts)
ivstops6 <- ivreg(log(Stops+1) ~ log(Placebo) + log(Ozone) + log(CO) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Agency) | . - log(Placebo) + log(PM_PY) + XXX, data = counts)
ivstops7 <- ivreg(log(Stops+1) ~ log(Placebo) + log(Ozone) + log(CO) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(Jurisdiction) | . - log(Placebo) + log(PM_PY) + XXX, data = counts)
ivstops8 <- ivreg(log(Stops+1) ~ log(Placebo) + log(Ozone) + log(CO) + TEMP + I(TEMP^2) + I(TEMP^3) + PRCP + factor(WEEK) + factor(DOW) + factor(Agency) | . - log(Placebo) + log(PM_PY) + XXX, data = counts)

ivstops1x <- coeftest(ivstops1, vcov = vcovCL(ivstops1, type = 'HC1'))
ivstops2x <- coeftest(ivstops2, vcov = vcovCL(ivstops2, type = 'HC1'))
ivstops3x <- coeftest(ivstops3, vcov = vcovCL(ivstops3, type = 'HC1'))
ivstops4x <- coeftest(ivstops4, vcov = vcovCL(ivstops4, type = 'HC1'))
ivstops5x <- coeftest(ivstops5, vcov = vcovCL(ivstops5, type = 'HC1'))
ivstops6x <- coeftest(ivstops6, vcov = vcovCL(ivstops6, type = 'HC1'))
ivstops7x <- coeftest(ivstops7, vcov = vcovCL(ivstops7, type = 'HC1'))
ivstops8x <- coeftest(ivstops8, vcov = vcovCL(ivstops8, type = 'HC1'))

write.csv(stargazer(ivstops1, ivstops2, ivstops3, ivstops4, ivstops5, ivstops6, ivstops7, ivstops8, type = 'text',
                    omit = c('WEEK', 'DOW', 'Agency', 'Jurisdiction')), paste0(direc, 'results/pm_iv_counts_unadj_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(ivstops1x, ivstops2x, ivstops3x, ivstops4x, ivstops5x, ivstops6x, ivstops7x, ivstops8x, type = 'text',
                    omit = c('WEEK', 'DOW', 'Agency', 'Jurisdiction')), paste0(direc, 'results/pm_iv_counts_placebo_pm.txt'), row.names = FALSE)

sum_stops1 <- summary(ivstops1, diagnostics = TRUE)
sum_stops2 <- summary(ivstops2, diagnostics = TRUE)
sum_stops3 <- summary(ivstops3, diagnostics = TRUE)
sum_stops4 <- summary(ivstops4, diagnostics = TRUE)
sum_stops5 <- summary(ivstops5, diagnostics = TRUE)

sum_stops6 <- summary(ivstops6, diagnostics = TRUE)
sum_stops7 <- summary(ivstops7, diagnostics = TRUE)
sum_stops8 <- summary(ivstops8, diagnostics = TRUE)

f_stats <- as.data.frame(cbind(c(sum_stops1[[12]][7], sum_stops2[[12]][7], sum_stops3[[12]][7], sum_stops4[[12]][7],
                                 sum_stops5[[12]][7], sum_stops6[[12]][7], sum_stops7[[12]][7], sum_stops8[[12]][7])))

colnames(f_stats) <- c('PM')

write.csv(f_stats, paste0(direc, 'results/f_stats_stops_placebo_pm.txt'), row.names = FALSE)

# Crashes regressions with lags

crash.new.ids <- sample.int(nrow(cdf), nrow(cdf))

cdf$crash.countP <- cdf$crash.count[crash.new.ids]
cdf$school.countP <- cdf$school.count[crash.new.ids]
cdf$school2.countP <- cdf$school2.count[crash.new.ids]
cdf$work.countP <- cdf$work.count[crash.new.ids]
cdf$unrestrained.countP <- cdf$unrestrained.count[crash.new.ids]
cdf$alcohol.countP <- cdf$alcohol.count[crash.new.ids]
cdf$animal.countP <- cdf$animal.count[crash.new.ids]
cdf$pedestrian.countP <- cdf$pedestrian.count[crash.new.ids]
cdf$hitrun.countP <- cdf$hitrun.count[crash.new.ids]
cdf$night.countP <- cdf$night.count[crash.new.ids]
cdf$young.countP <- cdf$young.count[crash.new.ids]
cdf$senior.countP <- cdf$senior.count[crash.new.ids]

c1 <- lm(crash.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c2 <- lm(school.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c3 <- lm(school2.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c4 <- lm(work.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c5 <- lm(unrestrained.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c6 <- lm(alcohol.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c7 <- lm(animal.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c8 <- lm(pedestrian.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c9 <- lm(hitrun.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c10 <- lm(night.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c11 <- lm(young.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c12 <- lm(senior.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)

ic1 <- ivreg(crash.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic2 <- ivreg(school.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic3 <- ivreg(school2.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic4 <- ivreg(work.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic5 <- ivreg(unrestrained.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic6 <- ivreg(alcohol.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic7 <- ivreg(animal.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic8 <- ivreg(pedestrian.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic9 <- ivreg(hitrun.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic10 <- ivreg(night.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic11 <- ivreg(young.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic12 <- ivreg(senior.countP ~ log(PM) + log(PM_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)

c1x <- coeftest(c1, vcov = vcovCL(c1), type = 'HC1')
c2x <- coeftest(c2, vcov = vcovCL(c2), type = 'HC1')
c3x <- coeftest(c3, vcov = vcovCL(c3), type = 'HC1')
c4x <- coeftest(c4, vcov = vcovCL(c4), type = 'HC1')
c5x <- coeftest(c5, vcov = vcovCL(c5), type = 'HC1')
c6x <- coeftest(c6, vcov = vcovCL(c6), type = 'HC1')
c7x <- coeftest(c7, vcov = vcovCL(c7), type = 'HC1')
c8x <- coeftest(c8, vcov = vcovCL(c8), type = 'HC1')
c9x <- coeftest(c9, vcov = vcovCL(c9), type = 'HC1')
c10x <- coeftest(c10, vcov = vcovCL(c10), type = 'HC1')
c11x <- coeftest(c11, vcov = vcovCL(c11), type = 'HC1')
c12x <- coeftest(c12, vcov = vcovCL(c12), type = 'HC1')

ic1x <- coeftest(ic1, vcov = vcovCL(ic1), type = 'HC1')
ic2x <- coeftest(ic2, vcov = vcovCL(ic2), type = 'HC1')
ic3x <- coeftest(ic3, vcov = vcovCL(ic3), type = 'HC1')
ic4x <- coeftest(ic4, vcov = vcovCL(ic4), type = 'HC1')
ic5x <- coeftest(ic5, vcov = vcovCL(ic5), type = 'HC1')
ic6x <- coeftest(ic6, vcov = vcovCL(ic6), type = 'HC1')
ic7x <- coeftest(ic7, vcov = vcovCL(ic7), type = 'HC1')
ic8x <- coeftest(ic8, vcov = vcovCL(ic8), type = 'HC1')
ic9x <- coeftest(ic9, vcov = vcovCL(ic9), type = 'HC1')
ic10x <- coeftest(ic10, vcov = vcovCL(ic10), type = 'HC1')
ic11x <- coeftest(ic11, vcov = vcovCL(ic11), type = 'HC1')
ic12x <- coeftest(ic12, vcov = vcovCL(ic12), type = 'HC1')

# Saving results

write.csv(stargazer(ic1x, ic2x, ic3x, ic4x, ic5x, ic6x, ic7x, ic8x, ic9x, ic10x, ic11x, ic12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_crashes_lag_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(c1x, c2x, c3x, c4x, c5x, c6x, c7x, c8x, c9x, c10x, c11x, c12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_crashes_lag_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, ic9, ic10, ic11, ic12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_unadj_crashes_lag_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_crashes_lag_placebo_outcomes.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum1 <- summary(ic1, diagnostics = TRUE)
sum2 <- summary(ic2, diagnostics = TRUE)
sum3 <- summary(ic3, diagnostics = TRUE)
sum4 <- summary(ic4, diagnostics = TRUE)
sum5 <- summary(ic5, diagnostics = TRUE)
sum6 <- summary(ic6, diagnostics = TRUE)
sum7 <- summary(ic7, diagnostics = TRUE)
sum8 <- summary(ic8, diagnostics = TRUE)
sum9 <- summary(ic9, diagnostics = TRUE)
sum10 <- summary(ic10, diagnostics = TRUE)
sum11 <- summary(ic11, diagnostics = TRUE)
sum12 <- summary(ic12, diagnostics = TRUE)

f_stats <- as.data.frame(c(sum1[[12]][7], sum2[[12]][7], sum3[[12]][7], sum4[[12]][7], sum5[[12]][7], sum6[[12]][7],
                           sum7[[12]][7], sum8[[12]][7], sum9[[12]][7], sum10[[12]][7], sum11[[12]][7], sum12[[12]][7]))

colnames(f_stats) <- c('PM')

write.csv(f_stats, paste0(direc, 'results/f_stats_crashes_lag_placebo_outcomes.txt'), row.names = FALSE)

# Crashes regressions without lags

c1 <- lm(crash.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c2 <- lm(school.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c3 <- lm(school2.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c4 <- lm(work.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c5 <- lm(unrestrained.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c6 <- lm(alcohol.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c7 <- lm(animal.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c8 <- lm(pedestrian.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c9 <- lm(hitrun.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c10 <- lm(night.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c11 <- lm(young.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c12 <- lm(senior.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)

ic1 <- ivreg(crash.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic2 <- ivreg(school.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic3 <- ivreg(school2.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic4 <- ivreg(work.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic5 <- ivreg(unrestrained.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic6 <- ivreg(alcohol.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic7 <- ivreg(animal.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic8 <- ivreg(pedestrian.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic9 <- ivreg(hitrun.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic10 <- ivreg(night.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic11 <- ivreg(young.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)
ic12 <- ivreg(senior.countP ~ log(PM) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(PM) + log(PM_PY) + XXX, data = cdf)

c1x <- coeftest(c1, vcov = vcovCL(c1), type = 'HC1')
c2x <- coeftest(c2, vcov = vcovCL(c2), type = 'HC1')
c3x <- coeftest(c3, vcov = vcovCL(c3), type = 'HC1')
c4x <- coeftest(c4, vcov = vcovCL(c4), type = 'HC1')
c5x <- coeftest(c5, vcov = vcovCL(c5), type = 'HC1')
c6x <- coeftest(c6, vcov = vcovCL(c6), type = 'HC1')
c7x <- coeftest(c7, vcov = vcovCL(c7), type = 'HC1')
c8x <- coeftest(c8, vcov = vcovCL(c8), type = 'HC1')
c9x <- coeftest(c9, vcov = vcovCL(c9), type = 'HC1')
c10x <- coeftest(c10, vcov = vcovCL(c10), type = 'HC1')
c11x <- coeftest(c11, vcov = vcovCL(c11), type = 'HC1')
c12x <- coeftest(c12, vcov = vcovCL(c12), type = 'HC1')

ic1x <- coeftest(ic1, vcov = vcovCL(ic1), type = 'HC1')
ic2x <- coeftest(ic2, vcov = vcovCL(ic2), type = 'HC1')
ic3x <- coeftest(ic3, vcov = vcovCL(ic3), type = 'HC1')
ic4x <- coeftest(ic4, vcov = vcovCL(ic4), type = 'HC1')
ic5x <- coeftest(ic5, vcov = vcovCL(ic5), type = 'HC1')
ic6x <- coeftest(ic6, vcov = vcovCL(ic6), type = 'HC1')
ic7x <- coeftest(ic7, vcov = vcovCL(ic7), type = 'HC1')
ic8x <- coeftest(ic8, vcov = vcovCL(ic8), type = 'HC1')
ic9x <- coeftest(ic9, vcov = vcovCL(ic9), type = 'HC1')
ic10x <- coeftest(ic10, vcov = vcovCL(ic10), type = 'HC1')
ic11x <- coeftest(ic11, vcov = vcovCL(ic11), type = 'HC1')
ic12x <- coeftest(ic12, vcov = vcovCL(ic12), type = 'HC1')

# Saving results

write.csv(stargazer(ic1x, ic2x, ic3x, ic4x, ic5x, ic6x, ic7x, ic8x, ic9x, ic10x, ic11x, ic12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_crashes_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(c1x, c2x, c3x, c4x, c5x, c6x, c7x, c8x, c9x, c10x, c11x, c12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_crashes_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, ic9, ic10, ic11, ic12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_unadj_crashes_placebo_outcomes.txt'), row.names = FALSE)

write.csv(stargazer(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_crashes_placebo_outcomes.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum1 <- summary(ic1, diagnostics = TRUE)
sum2 <- summary(ic2, diagnostics = TRUE)
sum3 <- summary(ic3, diagnostics = TRUE)
sum4 <- summary(ic4, diagnostics = TRUE)
sum5 <- summary(ic5, diagnostics = TRUE)
sum6 <- summary(ic6, diagnostics = TRUE)
sum7 <- summary(ic7, diagnostics = TRUE)
sum8 <- summary(ic8, diagnostics = TRUE)
sum9 <- summary(ic9, diagnostics = TRUE)
sum10 <- summary(ic10, diagnostics = TRUE)
sum11 <- summary(ic11, diagnostics = TRUE)
sum12 <- summary(ic12, diagnostics = TRUE)

f_stats <- as.data.frame(c(sum1[[12]][7], sum2[[12]][7], sum3[[12]][7], sum4[[12]][7], sum5[[12]][7], sum6[[12]][7],
                           sum7[[12]][7], sum8[[12]][7], sum9[[12]][7], sum10[[12]][7], sum11[[12]][7], sum12[[12]][7]))

colnames(f_stats) <- c('PM')

write.csv(f_stats, paste0(direc, 'results/f_stats_crashes_placebo_outcomes.txt'), row.names = FALSE)

# Crashes regressions with lags

cdf$Placebo <- cdf$PM[crash.new.ids]
cdf$Placebo_Lag <- cdf$PM_PY[crash.new.ids]

c1 <- lm(crash.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c2 <- lm(school.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c3 <- lm(school2.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c4 <- lm(work.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c5 <- lm(unrestrained.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c6 <- lm(alcohol.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c7 <- lm(animal.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c8 <- lm(pedestrian.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c9 <- lm(hitrun.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c10 <- lm(night.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c11 <- lm(young.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c12 <- lm(senior.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)

ic1 <- ivreg(crash.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic2 <- ivreg(school.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic3 <- ivreg(school2.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic4 <- ivreg(work.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic5 <- ivreg(unrestrained.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic6 <- ivreg(alcohol.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic7 <- ivreg(animal.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic8 <- ivreg(pedestrian.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic9 <- ivreg(hitrun.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic10 <- ivreg(night.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic11 <- ivreg(young.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic12 <- ivreg(senior.count ~ log(Placebo) + log(Placebo_Lag) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)

c1x <- coeftest(c1, vcov = vcovCL(c1), type = 'HC1')
c2x <- coeftest(c2, vcov = vcovCL(c2), type = 'HC1')
c3x <- coeftest(c3, vcov = vcovCL(c3), type = 'HC1')
c4x <- coeftest(c4, vcov = vcovCL(c4), type = 'HC1')
c5x <- coeftest(c5, vcov = vcovCL(c5), type = 'HC1')
c6x <- coeftest(c6, vcov = vcovCL(c6), type = 'HC1')
c7x <- coeftest(c7, vcov = vcovCL(c7), type = 'HC1')
c8x <- coeftest(c8, vcov = vcovCL(c8), type = 'HC1')
c9x <- coeftest(c9, vcov = vcovCL(c9), type = 'HC1')
c10x <- coeftest(c10, vcov = vcovCL(c10), type = 'HC1')
c11x <- coeftest(c11, vcov = vcovCL(c11), type = 'HC1')
c12x <- coeftest(c12, vcov = vcovCL(c12), type = 'HC1')

ic1x <- coeftest(ic1, vcov = vcovCL(ic1), type = 'HC1')
ic2x <- coeftest(ic2, vcov = vcovCL(ic2), type = 'HC1')
ic3x <- coeftest(ic3, vcov = vcovCL(ic3), type = 'HC1')
ic4x <- coeftest(ic4, vcov = vcovCL(ic4), type = 'HC1')
ic5x <- coeftest(ic5, vcov = vcovCL(ic5), type = 'HC1')
ic6x <- coeftest(ic6, vcov = vcovCL(ic6), type = 'HC1')
ic7x <- coeftest(ic7, vcov = vcovCL(ic7), type = 'HC1')
ic8x <- coeftest(ic8, vcov = vcovCL(ic8), type = 'HC1')
ic9x <- coeftest(ic9, vcov = vcovCL(ic9), type = 'HC1')
ic10x <- coeftest(ic10, vcov = vcovCL(ic10), type = 'HC1')
ic11x <- coeftest(ic11, vcov = vcovCL(ic11), type = 'HC1')
ic12x <- coeftest(ic12, vcov = vcovCL(ic12), type = 'HC1')

# Saving results

write.csv(stargazer(ic1x, ic2x, ic3x, ic4x, ic5x, ic6x, ic7x, ic8x, ic9x, ic10x, ic11x, ic12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_crashes_lag_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(c1x, c2x, c3x, c4x, c5x, c6x, c7x, c8x, c9x, c10x, c11x, c12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_crashes_lag_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, ic9, ic10, ic11, ic12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_unadj_crashes_lag_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_crashes_lag_placebo_pm.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum1 <- summary(ic1, diagnostics = TRUE)
sum2 <- summary(ic2, diagnostics = TRUE)
sum3 <- summary(ic3, diagnostics = TRUE)
sum4 <- summary(ic4, diagnostics = TRUE)
sum5 <- summary(ic5, diagnostics = TRUE)
sum6 <- summary(ic6, diagnostics = TRUE)
sum7 <- summary(ic7, diagnostics = TRUE)
sum8 <- summary(ic8, diagnostics = TRUE)
sum9 <- summary(ic9, diagnostics = TRUE)
sum10 <- summary(ic10, diagnostics = TRUE)
sum11 <- summary(ic11, diagnostics = TRUE)
sum12 <- summary(ic12, diagnostics = TRUE)

f_stats <- as.data.frame(c(sum1[[12]][7], sum2[[12]][7], sum3[[12]][7], sum4[[12]][7], sum5[[12]][7], sum6[[12]][7],
                           sum7[[12]][7], sum8[[12]][7], sum9[[12]][7], sum10[[12]][7], sum11[[12]][7], sum12[[12]][7]))

colnames(f_stats) <- c('PM')

write.csv(f_stats, paste0(direc, 'results/f_stats_crashes_lag_placebo_pm.txt'), row.names = FALSE)

# Crashes regressions without lags

c1 <- lm(crash.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c2 <- lm(school.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c3 <- lm(school2.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c4 <- lm(work.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c5 <- lm(unrestrained.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c6 <- lm(alcohol.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c7 <- lm(animal.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c8 <- lm(pedestrian.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c9 <- lm(hitrun.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c10 <- lm(night.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c11 <- lm(young.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)
c12 <- lm(senior.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC), data = cdf)

ic1 <- ivreg(crash.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic2 <- ivreg(school.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic3 <- ivreg(school2.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic4 <- ivreg(work.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic5 <- ivreg(unrestrained.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic6 <- ivreg(alcohol.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic7 <- ivreg(animal.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic8 <- ivreg(pedestrian.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic9 <- ivreg(hitrun.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic10 <- ivreg(night.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic11 <- ivreg(young.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)
ic12 <- ivreg(senior.count ~ log(Placebo) + factor(WEEK) + factor(DOW) + factor(JURISDIC) | . - log(Placebo) + log(PM_PY) + XXX, data = cdf)

c1x <- coeftest(c1, vcov = vcovCL(c1), type = 'HC1')
c2x <- coeftest(c2, vcov = vcovCL(c2), type = 'HC1')
c3x <- coeftest(c3, vcov = vcovCL(c3), type = 'HC1')
c4x <- coeftest(c4, vcov = vcovCL(c4), type = 'HC1')
c5x <- coeftest(c5, vcov = vcovCL(c5), type = 'HC1')
c6x <- coeftest(c6, vcov = vcovCL(c6), type = 'HC1')
c7x <- coeftest(c7, vcov = vcovCL(c7), type = 'HC1')
c8x <- coeftest(c8, vcov = vcovCL(c8), type = 'HC1')
c9x <- coeftest(c9, vcov = vcovCL(c9), type = 'HC1')
c10x <- coeftest(c10, vcov = vcovCL(c10), type = 'HC1')
c11x <- coeftest(c11, vcov = vcovCL(c11), type = 'HC1')
c12x <- coeftest(c12, vcov = vcovCL(c12), type = 'HC1')

ic1x <- coeftest(ic1, vcov = vcovCL(ic1), type = 'HC1')
ic2x <- coeftest(ic2, vcov = vcovCL(ic2), type = 'HC1')
ic3x <- coeftest(ic3, vcov = vcovCL(ic3), type = 'HC1')
ic4x <- coeftest(ic4, vcov = vcovCL(ic4), type = 'HC1')
ic5x <- coeftest(ic5, vcov = vcovCL(ic5), type = 'HC1')
ic6x <- coeftest(ic6, vcov = vcovCL(ic6), type = 'HC1')
ic7x <- coeftest(ic7, vcov = vcovCL(ic7), type = 'HC1')
ic8x <- coeftest(ic8, vcov = vcovCL(ic8), type = 'HC1')
ic9x <- coeftest(ic9, vcov = vcovCL(ic9), type = 'HC1')
ic10x <- coeftest(ic10, vcov = vcovCL(ic10), type = 'HC1')
ic11x <- coeftest(ic11, vcov = vcovCL(ic11), type = 'HC1')
ic12x <- coeftest(ic12, vcov = vcovCL(ic12), type = 'HC1')

# Saving results

write.csv(stargazer(ic1x, ic2x, ic3x, ic4x, ic5x, ic6x, ic7x, ic8x, ic9x, ic10x, ic11x, ic12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_crashes_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(c1x, c2x, c3x, c4x, c5x, c6x, c7x, c8x, c9x, c10x, c11x, c12x,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_crashes_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, ic9, ic10, ic11, ic12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK')),
          paste0(direc, 'results/pm_iv_unadj_crashes_placebo_pm.txt'), row.names = FALSE)

write.csv(stargazer(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
                    type = 'text', omit = c('DOW', 'JURISDIC', 'WEEK'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj_crashes_placebo_pm.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum1 <- summary(ic1, diagnostics = TRUE)
sum2 <- summary(ic2, diagnostics = TRUE)
sum3 <- summary(ic3, diagnostics = TRUE)
sum4 <- summary(ic4, diagnostics = TRUE)
sum5 <- summary(ic5, diagnostics = TRUE)
sum6 <- summary(ic6, diagnostics = TRUE)
sum7 <- summary(ic7, diagnostics = TRUE)
sum8 <- summary(ic8, diagnostics = TRUE)
sum9 <- summary(ic9, diagnostics = TRUE)
sum10 <- summary(ic10, diagnostics = TRUE)
sum11 <- summary(ic11, diagnostics = TRUE)
sum12 <- summary(ic12, diagnostics = TRUE)

f_stats <- as.data.frame(c(sum1[[12]][7], sum2[[12]][7], sum3[[12]][7], sum4[[12]][7], sum5[[12]][7], sum6[[12]][7],
                           sum7[[12]][7], sum8[[12]][7], sum9[[12]][7], sum10[[12]][7], sum11[[12]][7], sum12[[12]][7]))

colnames(f_stats) <- c('PM')

write.csv(f_stats, paste0(direc, 'results/f_stats_crashes_placebo_pm.txt'), row.names = FALSE)

