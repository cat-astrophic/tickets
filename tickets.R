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

' RANDOM FOREST APPROACH PREFERED ?!?! '

' WHAT ABOUT THE USE OF DAILY LAGS IN THE MODELS ?!?! '

' CAN ALSO TRY OZONE NOT +1 '

' ARE RESULTS ROBUST TO DROPPING XXX AS AN INSTRUMENT ?? '

' INTERACT RACE x POLLUTION AND GENDER x POLLUTION '


# Project directory

direc <- 'D:/tickets/'

# Read in the data

data <- read.csv(paste0(direc, 'data/data_raw.csv'))
counts <- read.csv(paste0(direc, 'data/counts_data.csv'))

# Creating additional variables

data$VEHICLE_SEARCHED <- ifelse(data$VEHICLE.SEARCHED == 'YES', 1, 0)
data$PERSON_SEARCHED <- ifelse(data$PERSON.SEARCHED == 'YES', 1, 0)
data$SEARCHED <- pmax(data$VEHICLE_SEARCHED, data$PERSON_SEARCHED)
data$SEARCHED_BOTH <- (data$VEHICLE_SEARCHED == data$PERSON_SEARCHED) * data$SEARCHED
data$ARRESTED <- ifelse(data$ACTION.TAKEN == 'ARREST', 1, 0)
data$TICKETED <- ifelse(data$ACTION.TAKEN == 'CITATION/SUMMONS', 1, 0)
data$FORCE <- ifelse(data$FORCE.USED.BY.OFFICER == 'YES', 1, 0)

# Subsetting to remove observations where the explanatory variables are ambiguous

data <- data %>% filter(ENGLISH.SPEAKING %in% c('YES', 'NO'))
data <- data %>% filter(GENDER %in% c('FEMALE', 'MALE'))
data <- data %>% filter(ETHNICITY %in% c('NOT HISPANIC OR LATINO', 'HISPANIC OR LATINO'))
data <- data %>% filter(RACE %in% c('WHITE', 'BLACK OR AFRICAN AMERICAN', 'ASIAN OR NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 'AMERICAN INDIAN OR ALASKA NATIVE'))
data <- data %>% filter(RESIDENCY %in% c('RESIDENT OF CITY/COUNTY OF STOP', 'OUT OF STATE RESIDENT', 'OTHER VIRGINIA JURISDICTION RESIDENT'))

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

peaks <- read_sf('D:/tickets/data/USA_Mountain_Peaks/v107/gsummit.gdb')
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
  ggtitle('Daily Mean Particulate Matter Levels by Week') + 
  ylab('Micrograms per Cubic Meter') +
  xlab('Day') +
  geom_path(aes(x = Day, y = `Mean Particulate Matter`, group = 1, color = Color), size = 1, alpha = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = 136) +
  geom_vline(xintercept = 205) +
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

# Save this version of data to minimize the intense crying

write.csv(data, paste0(direc, 'data/final_data.csv'), row.names = FALSE)

# OLS REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket.xxx <- lm(TICKETED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                 + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + factor(RESIDENCY) + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

ticket.x10 <- lm(TICKETED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                 + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + factor(RESIDENCY) + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xticket.xxx <- coeftest(ticket.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket.x10 <- coeftest(ticket.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested.xxx <- lm(ARRESTED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                   + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + factor(RESIDENCY) + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

arrested.x10 <- lm(ARRESTED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                   + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + factor(RESIDENCY) + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xarrested.xxx <- coeftest(arrested.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested.x10 <- coeftest(arrested.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx.xxx <- lm(ARRESTED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                    + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                    + factor(RESIDENCY) + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

arrestedx.x10 <- lm(ARRESTED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                    + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                    + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                    + factor(RESIDENCY) + factor(AGENCY.NAME)
                    + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

xarrestedx.xxx <- coeftest(arrestedx.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx.x10 <- coeftest(arrestedx.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched.xxx <- lm(SEARCHED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                   + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + factor(RESIDENCY) + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

searched.x10 <- lm(SEARCHED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                   + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                   + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                   + factor(RESIDENCY) + factor(AGENCY.NAME)
                   + factor(WEEK) + factor(DOW), data = data)

xsearched.xxx <- coeftest(searched.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched.x10 <- coeftest(searched.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle.xxx <- lm(VEHICLE_SEARCHED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                  + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                  + factor(RESIDENCY) + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

vehicle.x10 <- lm(VEHICLE_SEARCHED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                  + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                  + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                  + factor(RESIDENCY) + factor(AGENCY.NAME)
                  + factor(WEEK) + factor(DOW), data = data)

xvehicle.xxx <- coeftest(vehicle.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle.x10 <- coeftest(vehicle.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver.xxx <- lm(PERSON_SEARCHED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                 + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + factor(RESIDENCY) + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

driver.x10 <- lm(PERSON_SEARCHED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                 + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + factor(RESIDENCY) + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data)

xdriver.xxx <- coeftest(driver.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver.x10 <- coeftest(driver.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force.xxx <- lm(FORCE ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                + factor(RESIDENCY) + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

force.x10 <- lm(FORCE ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                + factor(RESIDENCY) + factor(AGENCY.NAME)
                + factor(WEEK) + factor(DOW), data = data)

xforce.xxx <- coeftest(force.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce.x10 <- coeftest(force.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex.xxx <- lm(FORCE ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                 + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + factor(RESIDENCY) + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

forcex.x10 <- lm(FORCE ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                 + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                 + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                 + factor(RESIDENCY) + factor(AGENCY.NAME)
                 + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

xforcex.xxx <- coeftest(forcex.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex.x10 <- coeftest(forcex.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# IV REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(RESIDENCY) + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

ticket_iv.x10 <- ivreg(TICKETED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(RESIDENCY) + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                         + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(RESIDENCY) + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

arrested_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                         + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(RESIDENCY) + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                          + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(RESIDENCY) + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                          + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(RESIDENCY) + factor(AGENCY.NAME)
                          + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                         + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(RESIDENCY) + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

searched_iv.x10 <- ivreg(SEARCHED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                         + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(RESIDENCY) + factor(AGENCY.NAME)
                         + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                        + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(RESIDENCY) + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                        + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(RESIDENCY) + factor(AGENCY.NAME)
                        + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(RESIDENCY) + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(RESIDENCY) + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                      + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(RESIDENCY) + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data)

force_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                      + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(RESIDENCY) + factor(AGENCY.NAME)
                      + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ log(PM) + log(Ozone+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(RESIDENCY) + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM) + log(PM_PY) + XXX, data = data[which(data$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ log(PM10) + log(Ozone+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(RESIDENCY) + factor(AGENCY.NAME)
                       + factor(WEEK) + factor(DOW) | . - log(PM10) + log(PM10_PY) + XXX, data = data[which(data$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
          xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE')),
          paste0(direc, 'results/pm_iv.txt'), row.names = FALSE)

write.csv(stargazer(xticket.xxx, xarrested.xxx, xarrestedx.xxx, xsearched.xxx, xvehicle.xxx, xdriver.xxx, xforce.xxx, xforcex.xxx,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
          force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE')),
          paste0(direc, 'results/pm_iv_unadj.txt'), row.names = FALSE)

write.csv(stargazer(ticket.xxx, arrested.xxx, arrestedx.xxx, searched.xxx, vehicle.xxx, driver.xxx, force.xxx, forcex.xxx,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
          xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE')),
          paste0(direc, 'results/pm10_iv.txt'), row.names = FALSE)

write.csv(stargazer(xticket.x10, xarrested.x10, xarrestedx.x10, xsearched.x10, xvehicle.x10, xdriver.x10, xforce.x10, xforcex.x10,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
          force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE')),
          paste0(direc, 'results/pm10_iv_unadj.txt'), row.names = FALSE)

write.csv(stargazer(ticket.x10, arrested.x10, arrestedx.x10, searched.x10, vehicle.x10, driver.x10, force.x10, forcex.x10,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm10_ols_unadj.txt'), row.names = FALSE)

# Viewing results

stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
          xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'))

stargazer(xticket.xxx, xarrested.xxx, xarrestedx.xxx, xsearched.xxx, xvehicle.xxx, xdriver.xxx, xforce.xxx, xforcex.xxx,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f'))

stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
          force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'))

stargazer(ticket.xxx, arrested.xxx, arrestedx.xxx, searched.xxx, vehicle.xxx, driver.xxx, force.xxx, forcex.xxx,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f'))

stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
          xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'))

stargazer(xticket.x10, xarrested.x10, xarrestedx.x10, xsearched.x10, xvehicle.x10, xdriver.x10, xforce.x10, xforcex.x10,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f'))

stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
          force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'))

stargazer(ticket.x10, arrested.x10, arrestedx.x10, searched.x10, vehicle.x10, driver.x10, force.x10, forcex.x10,
          type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'), omit.stat = c('ser', 'f'))

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




















# Add weather data to counts

counts <- left_join(counts, data[,which(colnames(data) %in% c('STOP_DATE', 'JURISDICTION', 'Ozone', 'CO', 'TEMP', 'PRCP', 'PM_PY', 'PD', 'WEEK', 'DAY', 'Ozone_PY', 'CO_PY'))], by = c('DATE2' = 'STOP_DATE', 'Jurisdiction' = 'JURISDICTION'))
counts <- counts[!duplicated(counts),]

# Counts analyses

# Stops

stops1 <- lm(log(Stops+1) ~ log(PM+1) + TEMP + PRCP, data = counts[which(counts$DAY %in% c(136:205)),])
stops2 <- lm(log(Stops+1) ~ log(PM+1) + TEMP + PRCP + factor(Agency), data = counts[which(counts$DAY %in% c(136:205)),])

ivstops1 <- ivreg(log(Stops+1) ~ log(PM+1) + TEMP + PRCP | . - log(PM+1) + log(PD+1), data = counts[which(counts$DAY %in% c(136:205)),])
ivstops2 <- ivreg(log(Stops+1) ~ log(PM+1) + TEMP + PRCP + factor(Agency) | . - log(PM+1) + log(PMD+1), data = counts[which(counts$DAY %in% c(136:205)),])

stops1x <- coeftest(stops1, vcov = vcovCL(stops1, type = 'HC0'))
stops2x <- coeftest(stops2, vcov = vcovCL(stops2, type = 'HC0'))

ivstops1x <- coeftest(ivstops1, vcov = vcovCL(ivstops1, type = 'HC0'))
ivstops2x <- coeftest(ivstops2, vcov = vcovCL(ivstops2, type = 'HC0'))

stargazer(stops1, ivstops1, stops2, ivstops2, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

stargazer(stops1x, ivstops1x, stops2x, ivstops2x, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

shit.ass <- summary(ivstops2, diagnostics = TRUE)






























# Tickets

tix1 <- lm(log(Tickets+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1), data = counts)
tix2 <- lm(log(Tickets+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency), data = counts)

ivtix1 <- ivreg(log(Tickets+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) | . - PM + IV_PM, data = counts)
ivtix2 <- ivreg(log(Tickets+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency) | . - PM + IV_PM, data = counts)

tix1x <- coeftest(tix1, vcov = vcovCL(tix1, type = 'HC0'))
tix2x <- coeftest(tix2, vcov = vcovCL(tix2, type = 'HC0'))

ivtix1x <- coeftest(ivtix1, vcov = vcovCL(ivtix1, type = 'HC0'))
ivtix2x <- coeftest(ivtix2, vcov = vcovCL(ivtix2, type = 'HC0'))

stargazer(tix1, ivtix1, tix2, ivtix2, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

stargazer(tix1x, ivtix1x, tix2x, ivtix2x, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

# Searches

sear1 <- lm(log(Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1), data = counts)
sear2 <- lm(log(Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency), data = counts)

ivsear1 <- ivreg(log(Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) | . - PM + IV_PM, data = counts)
ivsear2 <- ivreg(log(Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency) | . - PM + IV_PM, data = counts)

sear1x <- coeftest(sear1, vcov = vcovCL(sear1, type = 'HC0'))
sear2x <- coeftest(sear2, vcov = vcovCL(sear2, type = 'HC0'))

ivsear1x <- coeftest(ivsear1, vcov = vcovCL(ivsear1, type = 'HC0'))
ivsear2x <- coeftest(ivsear2, vcov = vcovCL(ivsear2, type = 'HC0'))

stargazer(sear1, ivsear1, sear2, ivsear2, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

stargazer(sear1x, ivsear1x, sear2x, ivsear2x, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

# Vehicle Searches

vsear1 <- lm(log(Vehicles_Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1), data = counts)
vsear2 <- lm(log(Vehicles_Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency), data = counts)

ivvsear1 <- ivreg(log(Vehicles_Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) | . - PM + IV_PM, data = counts)
ivvsear2 <- ivreg(log(Vehicles_Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency) | . - PM + IV_PM, data = counts)

vsear1x <- coeftest(vsear1, vcov = vcovCL(vsear1, type = 'HC0'))
vsear2x <- coeftest(vsear2, vcov = vcovCL(vsear2, type = 'HC0'))

ivvsear1x <- coeftest(ivvsear1, vcov = vcovCL(ivvsear1, type = 'HC0'))
ivvsear2x <- coeftest(ivvsear2, vcov = vcovCL(ivvsear2, type = 'HC0'))

stargazer(vsear1, ivvsear1, vsear2, ivvsear2, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

stargazer(vsear1x, ivvsear1x, vsear2x, ivvsear2x, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

# Person Searches

psear1 <- lm(log(Persons_Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1), data = counts)
psear2 <- lm(log(Persons_Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency), data = counts)

ivpsear1 <- ivreg(log(Persons_Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) | . - PM + IV_PM, data = counts)
ivpsear2 <- ivreg(log(Persons_Searched+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency) | . - PM + IV_PM, data = counts)

psear1x <- coeftest(psear1, vcov = vcovCL(psear1, type = 'HC0'))
psear2x <- coeftest(psear2, vcov = vcovCL(psear2, type = 'HC0'))

ivpsear1x <- coeftest(ivpsear1, vcov = vcovCL(ivpsear1, type = 'HC0'))
ivpsear2x <- coeftest(ivpsear2, vcov = vcovCL(ivpsear2, type = 'HC0'))

stargazer(psear1, ivpsear1, psear2, ivpsear2, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

stargazer(psear1x, ivpsear1x, psear2x, ivpsear2x, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

# Arrests

arrest1 <- lm(log(Arrests+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1), data = counts)
arrest2 <- lm(log(Arrests+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency), data = counts)

ivarrest1 <- ivreg(log(Arrests+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) | . - PM + IV_PM, data = counts)
ivarrest2 <- ivreg(log(Arrests+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency) | . - PM + IV_PM, data = counts)

arrest1x <- coeftest(arrest1, vcov = vcovCL(arrest1, type = 'HC0'))
arrest2x <- coeftest(arrest2, vcov = vcovCL(arrest2, type = 'HC0'))

ivarrest1x <- coeftest(ivarrest1, vcov = vcovCL(ivarrest1, type = 'HC0'))
ivarrest2x <- coeftest(ivarrest2, vcov = vcovCL(ivarrest2, type = 'HC0'))

stargazer(arrest1, ivarrest1, arrest2, ivarrest2, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

stargazer(arrest1x, ivarrest1x, arrest2x, ivarrest2x, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

# Force Used

force1 <- lm(log(Force_Used+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1), data = counts)
force2 <- lm(log(Force_Used+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency), data = counts)

ivforce1 <- ivreg(log(Force_Used+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) | . - PM + IV_PM, data = counts)
ivforce2 <- ivreg(log(Force_Used+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + factor(Agency) | . - PM + IV_PM, data = counts)

force1x <- coeftest(force1, vcov = vcovCL(force1, type = 'HC0'))
force2x <- coeftest(force2, vcov = vcovCL(force2, type = 'HC0'))

ivforce1x <- coeftest(ivforce1, vcov = vcovCL(ivforce1, type = 'HC0'))
ivforce2x <- coeftest(ivforce2, vcov = vcovCL(ivforce2, type = 'HC0'))

stargazer(force1, ivforce1, force2, ivforce2, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

stargazer(force1x, ivforce1x, force2x, ivforce2x, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

# Force Used | Arrested

aforce1 <- lm(log(Force_Used_Arrested+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + log(Arrests+1), data = counts)
aforce2 <- lm(log(Force_Used_Arrested+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + log(Arrests+1) + factor(Agency), data = counts)

ivaforce1 <- ivreg(log(Force_Used_Arrested+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + log(Arrests+1) | . - PM + IV_PM, data = counts)
ivaforce2 <- ivreg(log(Force_Used_Arrested+1) ~ PM + log(Ozone+1) + log(CO+1) + TEMP + PRCP + log(Stops+1) + log(Arrests+1) + factor(Agency) | . - PM + IV_PM, data = counts)

aforce1x <- coeftest(aforce1, vcov = vcovCL(aforce1, type = 'HC0'))
aforce2x <- coeftest(aforce2, vcov = vcovCL(aforce2, type = 'HC0'))

ivaforce1x <- coeftest(ivaforce1, vcov = vcovCL(ivaforce1, type = 'HC0'))
ivaforce2x <- coeftest(ivaforce2, vcov = vcovCL(ivaforce2, type = 'HC0'))

stargazer(aforce1, ivaforce1, aforce2, ivaforce2, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

stargazer(aforce1x, ivaforce1x, aforce2x, ivaforce2x, omit = c('Agency', 'WEEK', 'Jurisdiction'), type = 'text')

# Save results

write.csv(stargazer(stops2, tix2, sear2, vsear2, psear2, arrest2, force2, aforce2, omit = c('Agency', 'WEEK', 'Jurisdiction'),
                    omit.stat = c('f', 'ser')), paste0(direc, 'results/counts_ols_unadj.txt'), row.names = FALSE)

write.csv(stargazer(stops2x, tix2x, sear2x, vsear2x, psear2x, arrest2x, force2x, aforce2x, omit = c('Agency', 'WEEK', 'Jurisdiction'),
                    omit.stat = c('f', 'ser')), paste0(direc, 'results/counts_ols_adj.txt'), row.names = FALSE)

write.csv(stargazer(ivstops2, ivtix2, ivsear2, ivvsear2, ivpsear2, ivarrest2, ivforce2, ivaforce2, omit = c('Agency', 'WEEK', 'Jurisdiction'),
                    omit.stat = c('f', 'ser')), paste0(direc, 'results/counts_iv_unadj.txt'), row.names = FALSE)

write.csv(stargazer(ivstops2x, ivtix2x, ivsear2x, ivvsear2x, ivpsear2x, ivarrest2x, ivforce2x, ivaforce2x, omit = c('Agency', 'WEEK', 'Jurisdiction'),
                    omit.stat = c('f', 'ser')), paste0(direc, 'results/counts_iv_adj.txt'), row.names = FALSE)

# Checking results together

stargazer(stops1, tix1, sear1, vsear1, psear1, arrest1, force1, aforce1, type = 'text', omit = c('Agency', 'WEEK', 'Jurisdiction'), omit.stat = c('f', 'ser'))

stargazer(stops1x, tix1x, sear1x, vsear1x, psear1x, arrest1x, force1x, aforce1x, type = 'text', omit = c('Agency', 'WEEK', 'Jurisdiction'), omit.stat = c('f', 'ser'))

stargazer(ivstops1, ivtix1, ivsear1, ivvsear1, ivpsear1, ivarrest1, ivforce1, ivaforce1, type = 'text', omit = c('Agency', 'WEEK', 'Jurisdiction'), omit.stat = c('f', 'ser'))

stargazer(ivstops1x, ivtix1x, ivsear1x, ivvsear1x, ivpsear1x, ivarrest1x, ivforce1x, ivaforce1x, type = 'text', omit = c('Agency', 'WEEK', 'Jurisdiction'), omit.stat = c('f', 'ser'))

stargazer(stops2, tix2, sear2, vsear2, psear2, arrest2, force2, aforce2, type = 'text', omit = c('Agency', 'WEEK', 'Jurisdiction'), omit.stat = c('f', 'ser'))

stargazer(stops2x, tix2x, sear2x, vsear2x, psear2x, arrest2x, force2x, aforce2x, type = 'text', omit = c('Agency', 'WEEK', 'Jurisdiction'), omit.stat = c('f', 'ser'))

stargazer(ivstops2, ivtix2, ivsear2, ivvsear2, ivpsear2, ivarrest2, ivforce2, ivaforce2, type = 'text', omit = c('Agency', 'WEEK', 'Jurisdiction'), omit.stat = c('f', 'ser'))

stargazer(ivstops2x, ivtix2x, ivsear2x, ivvsear2x, ivpsear2x, ivarrest2x, ivforce2x, ivaforce2x, type = 'text', omit = c('Agency', 'WEEK', 'Jurisdiction'), omit.stat = c('f', 'ser'))

# First stage F-statistics from above

sum_stops <- summary(ivstops1, diagnostics = TRUE)
sum_tix <- summary(ivtix1, diagnostics = TRUE)
sum_sear <- summary(ivsear1, diagnostics = TRUE)
sum_vsear <- summary(ivvsear1, diagnostics = TRUE)
sum_psear <- summary(ivpsear1, diagnostics = TRUE)
sum_arrest <- summary(ivarrest1, diagnostics = TRUE)
sum_force <- summary(ivforce1, diagnostics = TRUE)
sum_aforce <- summary(ivaforce1, diagnostics = TRUE)

sum_stops2 <- summary(ivstops2, diagnostics = TRUE)
sum_tix2 <- summary(ivtix2, diagnostics = TRUE)
sum_sear2 <- summary(ivsear2, diagnostics = TRUE)
sum_vsear2 <- summary(ivvsear2, diagnostics = TRUE)
sum_psear2 <- summary(ivpsear2, diagnostics = TRUE)
sum_arrest2 <- summary(ivarrest2, diagnostics = TRUE)
sum_force2 <- summary(ivforce2, diagnostics = TRUE)
sum_aforce2 <- summary(ivaforce2, diagnostics = TRUE)

f_stats2 <- as.data.frame(cbind(c(sum_stops[[12]][7], sum_tix[[12]][7], sum_sear[[12]][7], sum_vsear[[12]][7],
                                  sum_psear[[12]][7], sum_arrest[[12]][7], sum_force[[12]][7], sum_aforce[[12]][7]),
                                c(sum_stops2[[12]][7], sum_tix2[[12]][7], sum_sear2[[12]][7], sum_vsear2[[12]][7],
                                  sum_psear2[[12]][7], sum_arrest2[[12]][7], sum_force2[[12]][7], sum_aforce2[[12]][7])))

colnames(f_stats2) <- c('Baseline', 'Main')

write.csv(f_stats2, paste0(direc, 'results/f_stats2.txt'), row.names = FALSE)

# Creating summary statistics tables

data$Female <- ifelse(data$GENDER == 'FEMALE', 1, 0)
data$Male <- ifelse(data$GENDER == 'FEMALE', 0, 1)
data$Hispanic <- ifelse(data$ETHNICITY == 'HISPANIC OR LATINO', 1, 0)
data$English <- ifelse(data$ENGLISH.SPEAKING == 'YES', 1, 0)
data$Asian <- ifelse(data$RACE == 'ASIAN OR NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 1, 0)
data$Black <- ifelse(data$RACE == 'BLACK OR AFRICAN AMERICAN', 1, 0)
data$Native <- ifelse(data$RACE == 'AMERICAN INDIAN OR ALASKA NATIVE', 1, 0)
data$White <- ifelse(data$RACE == 'WHITE', 1, 0)

dfsum <- data[,c(63, 62, 60, 58, 59, 64, 51, 52, 53, 55, 27, 45, 8, 73, 86, 90, 91, 92, 93, 88, 89)]

csum <- counts[,c(4:5, 15, 6:10, 12, 18:19, 16:17)]

colnames(dfsum) <- c('Ticketed', 'Arrested', 'Searched', 'Vehicle Searched',
                     'Person Searched', 'Force Used', 'Particulate Matter (PM2.5)',
                     'Particulate Matter (PM10)', 'Ozone', 'Carbon Monoxide',
                     'Temperature', 'Precipitation', 'Age', 'Minor', 'Female',
                     'Asian', 'Black', 'Native American', 'White', 'Hispanic', 'Speaks English')

colnames(csum) <- c('Stops', 'Tickets', 'Searches', 'Vehicle Searches', 'Person Searches',
                    'Arrests', 'Force Used', 'Force Used | Arrest Made', 'Particulate Matter (PM2.5)',
                    'Ozone', 'Carbon Monoxide', 'Temperature', 'Precipitaton')

datasummary_skim(dfsum, fmt = '%.3f')

datasummary_skim(csum, fmt = '%.3f')

# Re-run the IV analysis for same-county residents only

samesies <- data %>% filter(RESIDENCY == 'RESIDENT OF CITY/COUNTY OF STOP')

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED ~ PM + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM + IV_PM, data = samesies)

ticket_iv.x10 <- ivreg(TICKETED ~ PM10 + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM10 + IV_PM10, data = samesies)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xticket_iv.x10 <- coeftest(ticket_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ PM + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                         + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM + IV_PM, data = samesies)

arrested_iv.x10 <- ivreg(ARRESTED ~ PM10 + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                         + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM10 + IV_PM10, data = samesies)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrested_iv.x10 <- coeftest(arrested_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ PM + log(Ozone+1) + log(CO+1) + log(Stops+1) + log(Searched+1) + MINOR
                          + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(AGENCY.NAME)
                          + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM + IV_PM, data = samesies[which(samesies$SEARCHED == 1),])

arrestedx_iv.x10 <- ivreg(ARRESTED ~ PM10 + log(Ozone+1) + log(CO+1) + log(Stops+1) + log(Searched+1) + MINOR
                          + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                          + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                          + factor(AGENCY.NAME)
                          + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM10 + IV_PM10, data = samesies[which(samesies$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xarrestedx_iv.x10 <- coeftest(arrestedx_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ PM + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                         + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM + IV_PM, data = samesies)

searched_iv.x10 <- ivreg(SEARCHED ~ PM10 + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                         + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                         + factor(AGENCY.NAME)
                         + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM10 + IV_PM10, data = samesies)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xsearched_iv.x10 <- coeftest(searched_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ PM + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                        + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(AGENCY.NAME)
                        + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM + IV_PM, data = samesies)

vehicle_iv.x10 <- ivreg(VEHICLE_SEARCHED ~ PM10 + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                        + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                        + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                        + factor(AGENCY.NAME)
                        + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM10 + IV_PM10, data = samesies)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xvehicle_iv.x10 <- coeftest(vehicle_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ PM + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM + IV_PM, data = samesies)

driver_iv.x10 <- ivreg(PERSON_SEARCHED ~ PM10 + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                       + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM10 + IV_PM10, data = samesies)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xdriver_iv.x10 <- coeftest(driver_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ PM + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                      + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(AGENCY.NAME)
                      + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM + IV_PM, data = samesies)

force_iv.x10 <- ivreg(FORCE ~ PM10 + log(Ozone+1) + log(CO+1) + log(Stops+1) + MINOR
                      + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                      + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                      + factor(AGENCY.NAME)
                      + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM10 + IV_PM10, data = samesies)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforce_iv.x10 <- coeftest(force_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ PM + log(Ozone+1) + log(CO+1) + log(Stops+1) + log(Arrests+1) + MINOR
                       + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM + IV_PM, data = samesies[which(samesies$ARRESTED == 1),])

forcex_iv.x10 <- ivreg(FORCE ~ PM10 + log(Ozone+1) + log(CO+1) + log(Stops+1) + log(Arrests+1) + MINOR
                       + TEMP + PRCP + relevel(RACE, ref = 'WHITE') + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(REASON.FOR.STOP)
                       + factor(AGENCY.NAME)
                       + factor(JURISDICTION) + factor(WEEK) + factor(DOW) | . - PM10 + IV_PM10, data = samesies[which(samesies$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)
xforcex_iv.x10 <- coeftest(forcex_iv.x10, vcov = vcovCL, cluster = ~AGENCY.NAME)














# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE')),
          paste0(direc, 'results/pm_iv_local.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE')),
          paste0(direc, 'results/pm_iv_unadj_local.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
                    xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE')),
          paste0(direc, 'results/pm10_iv_local.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
                    force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE')),
          paste0(direc, 'results/pm10_iv_unadj_local.txt'), row.names = FALSE)

# Viewing results

stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
          xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'))

stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
          force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'))

stargazer(xticket_iv.x10, xarrested_iv.x10, xarrestedx_iv.x10, xsearched_iv.x10, xvehicle_iv.x10, xdriver_iv.x10,
          xforce_iv.x10, xforcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'))

stargazer(ticket_iv.x10, arrested_iv.x10, arrestedx_iv.x10, searched_iv.x10, vehicle_iv.x10, driver_iv.x10,
          force_iv.x10, forcex_iv.x10, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'WEEK', 'VIRGINIA.CRIME.CODE'))

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

write.csv(f_stats3, paste0(direc, 'results/f_stats_samesies.txt'), row.names = FALSE)





' Two instruments - one is same day of year (e.g., 5th Monday) and one is same date '
' in both cases the instrument includes the extra mountain peak and time based component '












' increased stops could be due to drivers OR due to cops being more strict with speeding and then changing their mind (an act against boredom?) '


# FIGURES

# Make figures of raw and conditional outcomes of each outcome by week and overlay the Canadian forest fires










wee <- c()
loc <- c()
dpm <- c()
dPM10 <- c()
dco <- c()
do3 <- c()
dno2 <- c()

for (i in 1:53) {
  
  tmpx <- data %>% filter(WEEK == i)
  
  for (you.knee in unique(data$JURISDICTION)) {
    
    tmp <- tmpx %>% filter(JURISDICTION == you.knee)
    wee <- c(wee, i)
    loc <- c(loc, which(unique(data$JURISDICTION) == you.knee))
    dpm <- c(dpm, mean(tmp$PM, na.rm = TRUE))
    dPM10 <- c(dPM10, mean(tmp$PM10, na.rm = TRUE))
    dco <- c(dco, mean(tmp$CO, na.rm = TRUE))
    do3 <- c(do3, mean(tmp$Ozone, na.rm = TRUE))
    dno2 <- c(dno2, mean(tmp$NO2, na.rm = TRUE)) 
    
  }
  
}

plot_df <- as.data.frame(cbind(wee, loc, dpm, dPM10, dco, do3, dno2))

plot_pm <- plot_df %>% filter(dpm > 0)
plot_PM10 <- plot_df %>% filter(dPM10 > 0)
plot_co <- plot_df %>% filter(dco > 0)
plot_o3 <- plot_df %>% filter(do3 > 0)
plot_no2 <- plot_df %>% filter(dno2 > 0)

for (l in unique(plot_pm$loc)) {
  
  poop <- plot_pm %>% filter(loc == l)
  
  if (l == unique(plot_pm$loc)[1]) {
    
    plot(poop$dpm, type = 'l', ylim = c(0, 60))
    
  } else {
    
    lines(poop$dpm)
    
  }
  
}

for (l in unique(plot_PM10$loc)) {
  
  poop <- plot_PM10 %>% filter(loc == l)
  
  if (l == unique(plot_PM10$loc)[1]) {
    
    plot(poop$dPM10, type = 'l', ylim = c(0, 75))
    
  } else {
    
    lines(poop$dPM10)
    
  }
  
}

for (l in unique(plot_co$loc)) {
  
  poop <- plot_co %>% filter(loc == l)
  
  if (l == unique(plot_co$loc)[1]) {
    
    plot(poop$dco, type = 'l', ylim = c(0, 1))
    
  } else {
    
    lines(poop$dco)
    
  }
  
}

for (l in unique(plot_o3$loc)) {
  
  poop <- plot_o3 %>% filter(loc == l)
  
  if (l == unique(plot_o3$loc)[1]) {
    
    plot(poop$do3, type = 'l', ylim = c(0, 60))
    
  } else {
    
    lines(poop$do3)
    
  }
  
}

for (l in unique(plot_no2$loc)) {
  
  poop <- plot_no2 %>% filter(loc == l)
  
  if (l == unique(plot_no2$loc)[1]) {
    
    plot(poop$dno2, type = 'l', ylim = c(0, 60))
    
  } else {
    
    lines(poop$dno2)
    
  }
  
}

# Make a map of counties for PM and PM10 included in sample colored by annual stops and out of sample with opacity zero


' the treatment needs to coincide with the actual fucking fires '
' it does because we dont know exactly when VA was affected by fires, but this is the best estimate (of the era) '



