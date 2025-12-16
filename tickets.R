# This script performs the econometric analysis for pollution and tickets

# Loading libraries

library(modelsummary)
library(stargazer)
library(sandwich)
library(ggplot2)
library(leaflet)
library(lmtest)
library(tigris)
library(dplyr)
library(MASS)
library(AER)
library(sf)

# Project directory

direc <- 'D:/tickets/'

# Read in the data

data <- read.csv(paste0(direc, 'data/raw_data.csv'))
counts <- read.csv(paste0(direc, 'data/counts_data.csv'))
wind <- read.csv(paste0(direc, 'data/wind_data.csv'))
pm <- read.csv(paste0(direc, 'data/daily_PM.csv'))
ax <- read.csv(paste0(direc, 'data/CrashData_test_6478750435646127290.csv'))
va <- read_sf(paste0(direc, 'data/shape/va_counties.shp'))

# Creating additional variables

data$VEHICLE_SEARCHED <- ifelse(data$VEHICLE.SEARCHED == 'Y', 1, 0)
data$PERSON_SEARCHED <- ifelse(data$PERSON.SEARCHED == 'Y', 1, 0)
data$SEARCHED <- pmax(data$VEHICLE_SEARCHED, data$PERSON_SEARCHED)
data$SEARCHED_BOTH <- (data$VEHICLE_SEARCHED == data$PERSON_SEARCHED) * data$SEARCHED
data$ARRESTED <- ifelse(data$ACTION.TAKEN == 'Arrest', 1, 0)
data$TICKETED <- ifelse(data$ACTION.TAKEN == 'Citation/Summons', 1, 0)
data$FORCE <- ifelse(data$FORCE.USED.BY.OFFICER == 'Y', 1, 0)

# Subsetting to remove observations where the explanatory variables are ambiguous

data <- data %>% filter(ENGLISH.SPEAKING %in% c('Y', 'N'))
data <- data %>% filter(GENDER %in% c('Female', 'Male'))
data <- data %>% filter(ETHNICITY %in% c('Non-Hispanic or Latino', 'Hispanic or Latino'))
data <- data %>% filter(RACE %in% c('White', 'Black or African American', 'Asian or Native Hawaiian or Other Pacific Islander', 'American Indian or Alaska Native'))
data <- data %>% filter(RESIDENCY %in% c('Resident of city/county of stop', 'Out of State resident', 'Other Virginia jurisdiction resident'))

# Create a week of year variable

shit <- function(day) {
  
  ref <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  a <- gregexpr('/', day)[[1]][1]
  b <- gregexpr('/', day)[[1]][2]
  
  val <- as.integer(substr(day, a+1, b-1)) + sum(ref[1:(as.integer(substr(day, 1, a-1)))]) - 31
  
  return(ceiling(val/7))
  
}

weak <- c()

for (i in 1:nrow(data)) {
  
  print(i)
  weak <- c(weak, shit(data$STOP_DATE[i]))
  
}

data$WEEK <- weak

# Flagging minors in the data

data$MINOR <- as.integer(data$AGE < 18)

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

# Combined searches data

counts$Searched <- counts$Vehicles_Searched + counts$Persons_Searched

# Adding counts data to data

data <- left_join(data, counts, by = c('STOP_DATE' = 'DATE2', 'JURISDICTION' = 'Jurisdiction', 'AGENCY.NAME' = 'Agency'))

# Adding PM data to data

pm$Date2 <- paste0(substr(pm$Date, 1, 4), '-', substr(pm$Date, 5, 6), '-', substr(pm$Date, 7, 8))

data <- left_join(data, pm, by = c('FIPS.x' = 'FIPS', 'Date' = 'Date2'))

# Adding wind data to data (some FIPS have multiple stations so I take the mean values)

w.date <- c()
w.speed <- c()
w.direc <- c()
w.fips <- c()

w.date.list <- unique(wind$Date)
w.fips.list <- unique(wind$FIPS)

for (f in w.fips.list) {
  
  print(f)
  tmp <- wind[which(wind$FIPS == f),]
  
  for (d in w.date.list) {
    
    tmp2 <- tmp[which(tmp$Date == d),]
    w.date <- c(w.date, d)
    w.speed <- c(w.speed, mean(tmp2$Speed, na.rm = TRUE))
    w.direc <- c(w.direc, mean(tmp2$Direction, na.rm = TRUE))
    w.fips <- c(w.fips, f)
    
  }
  
}

wind2 <- as.data.frame(cbind(w.date, w.fips, w.speed, w.direc))
colnames(wind2) <- c('Date', 'FIPS', 'Speed', 'Direction')
wind2$FIPS <- as.integer(wind2$FIPS)
wind2$Speed <- as.integer(wind2$Speed)
wind2$Direction <- as.integer(wind2$Direction)

data <- left_join(data, wind2, by = c('FIPS.x' = 'FIPS', 'Date' = 'Date'))

# Creating a day of week variable

dates <- sort(unique(data$Date))
dow <- c()

for (i in 1:nrow(data)) {
  
  print(i)
  dow <- c(dow, (which(dates == data$Date[i]) - 2) %% 7)
  
}

data$DOW <- dow

# Binary rain data

data$Rain <- ifelse(data$PRCP > 0 & data$PRCP < 90, 1, 0)

# Creating binary wind direction categories

data$Wind <- floor(data$Direction/90)

# Identifying high PM days

cutoff <- quantile(data$PM, c(.75), na.rm = TRUE)
data$PMX <- as.integer(data$PM > cutoff)

# Month FE

mos <- c()

for (i in 1:nrow(data)) {
  
  print(i)
  mos <- c(mos, strsplit(data$STOP_DATE[i], '/')[[1]][1])
  
}

data$Month <- mos

# Save this version of data to minimize the intense crying

write.csv(data, paste0(direc, 'data/final_data.csv'), row.names = FALSE)

# Re-leveling fixed effects

data <- within(data, RACE <- relevel(factor(RACE), ref = 'White'))
data <- within(data, ETHNICITY <- relevel(factor(ETHNICITY), ref = 'Non-Hispanic or Latino'))
data <- within(data, GENDER <- relevel(factor(GENDER), ref = 'Female'))
data <- within(data, RESIDENCY <- relevel(factor(RESIDENCY), ref = 'Resident of city/county of stop'))
data <- within(data, REASON.FOR.STOP <- relevel(factor(REASON.FOR.STOP), ref = 'Traffic Violation'))

# OLS BASELINE REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket.xxx <- lm(TICKETED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                 + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                 + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                 + factor(JURISDICTION), data = data)

xticket.xxx <- coeftest(ticket.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested.xxx <- lm(ARRESTED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                   + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                   + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                   + factor(JURISDICTION), data = data)

xarrested.xxx <- coeftest(arrested.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx.xxx <- lm(ARRESTED ~ PMX + log(Stops) + log(Searched) + factor(REASON.FOR.STOP) + MINOR
                    + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                    + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                    + factor(JURISDICTION), data = data[which(data$SEARCHED == 1),])

xarrestedx.xxx <- coeftest(arrestedx.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched.xxx <- lm(SEARCHED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                   + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                   + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                   + factor(JURISDICTION), data = data)

xsearched.xxx <- coeftest(searched.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle.xxx <- lm(VEHICLE_SEARCHED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                  + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                  + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                  + factor(JURISDICTION), data = data)

xvehicle.xxx <- coeftest(vehicle.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver.xxx <- lm(PERSON_SEARCHED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                 + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                 + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                 + factor(JURISDICTION), data = data)

xdriver.xxx <- coeftest(driver.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force.xxx <- lm(FORCE ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                + factor(JURISDICTION), data = data)

xforce.xxx <- coeftest(force.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex.xxx <- lm(FORCE ~ PMX + log(Stops) + log(Arrests) + factor(REASON.FOR.STOP) + MINOR
                 + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                 + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                 + factor(JURISDICTION), data = data[which(data$ARRESTED == 1),])

xforcex.xxx <- coeftest(forcex.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# IV REGRESSIONS

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                       + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                       + factor(JURISDICTION) | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                         + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                         + factor(JURISDICTION) | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ PMX + log(Stops) + log(Searched) + factor(REASON.FOR.STOP) + MINOR
                          + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                          + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                          + factor(JURISDICTION) | . - PMX + factor(Wind) + Speed + factor(WEEK) + factor(DOW), data = data[which(data$SEARCHED == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                         + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                         + factor(JURISDICTION) | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                        + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                        + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                        + factor(JURISDICTION) | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                       + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                       + factor(JURISDICTION) | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                      + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                      + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                      + factor(JURISDICTION) | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ PMX + log(Stops) + log(Arrests) + factor(REASON.FOR.STOP) + MINOR
                       + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month)
                       + factor(JURISDICTION) | . - PMX + factor(Wind) + Speed + factor(WEEK) + factor(DOW), data = data[which(data$ARRESTED == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'Month')),
          paste0(direc, 'results/pm_iv.txt'), row.names = FALSE)

write.csv(stargazer(xticket.xxx, xarrested.xxx, xarrestedx.xxx, xsearched.xxx, xvehicle.xxx, xdriver.xxx, xforce.xxx, xforcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'Month'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, omit = c('AGENCY.NAME', 'JURISDICTION', 'Month')),
          paste0(direc, 'results/pm_iv_tex.txt'), row.names = FALSE)

write.csv(stargazer(xticket.xxx, xarrested.xxx, xarrestedx.xxx, xsearched.xxx, xvehicle.xxx, xdriver.xxx, xforce.xxx, xforcex.xxx,
                    omit = c('AGENCY.NAME', 'JURISDICTION', 'Month'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_tex.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'Month')),
          paste0(direc, 'results/pm_iv_unadj.txt'), row.names = FALSE)

write.csv(stargazer(ticket.xxx, arrested.xxx, arrestedx.xxx, searched.xxx, vehicle.xxx, driver.xxx, force.xxx, forcex.xxx,
                    type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'Month'), omit.stat = c('ser', 'f')),
          paste0(direc, 'results/pm_ols_unadj.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_tix <- summary(ticket_iv.xxx, diagnostics = TRUE)
sum_arr <- summary(arrested_iv.xxx, diagnostics = TRUE)
sum_arx <- summary(arrestedx_iv.xxx, diagnostics = TRUE)
sum_sea <- summary(searched_iv.xxx, diagnostics = TRUE)
sum_veh <- summary(vehicle_iv.xxx, diagnostics = TRUE)
sum_dri <- summary(driver_iv.xxx, diagnostics = TRUE)
sum_for <- summary(force_iv.xxx, diagnostics = TRUE)
sum_fox <- summary(forcex_iv.xxx, diagnostics = TRUE)

f_stats <- as.data.frame(c(sum_tix[[12]][7], sum_arr[[12]][7], sum_arx[[12]][7], sum_sea[[12]][7],
                           sum_veh[[12]][7], sum_dri[[12]][7], sum_for[[12]][7], sum_fox[[12]][7]))

colnames(f_stats) <- c('F_STAT')

write.csv(f_stats, paste0(direc, 'results/f_stats.txt'), row.names = FALSE)

# Cleaning the crashes data

# Filter for 2022 data

ax <- ax %>% filter(Crash.Year == 2022)

# Create a date column that matches STOP_DATE

ax.dates <- c()

for (i in 1:nrow(ax)) {
  
  print(i)
  ax.dates <- c(ax.dates, strsplit(ax$Crash.Date[i], ' ')[[1]][1])
  
}

ax$Date <- ax.dates

# Convert ax to a spatial data.frame

ax <- ax %>% filter(is.na(x) == FALSE) %>% filter(is.na(y) == FALSE)
axe <- st_as_sf(ax, coords = c('x', 'y'))

# Linking observations to FIPS via tigris shapefiles

cos <- counties(state = 'VA')
cos$GEOID <- as.integer(cos$GEOID)

axe <- st_set_crs(axe, st_crs(cos))

inside <- st_within(axe$geometry, cos$geometry)

innards <- c()

for (i in 1:nrow(axe)) {
  
  print(i)
  innards <- c(innards, inside[[i]][1])
  
}

a.fips <- c()

for (i in 1:nrow(axe)) {
  
  print(i)
  a.fips <- c(a.fips, cos$GEOID[innards[i]])
  
}

axe$FIPS <- a.fips

# Bringing in PM and other data from data to axe

a.pm <- c()
a.wind <- c()
a.speed <- c()
a.temp <- c()
a.rain <- c()
a.week <- c()
a.dow <- c()

for (i in 1:nrow(axe)) {
  
  print(i)
  
  tmp <- data %>% filter(STOP_DATE == axe$Date[i])
  tmp <- tmp %>% filter(FIPS.x == axe$FIPS[i])
  
  a.pm <- c(a.pm, tmp$PMX[1])
  a.wind <- c(a.wind, tmp$Wind[1])
  a.speed <- c(a.speed, tmp$Speed[1])
  a.temp <- c(a.temp, tmp$TEMP[1])
  a.rain <- c(a.rain, tmp$Rain[1])
  a.week <- c(a.week, tmp$WEEK[1])
  a.dow <- c(a.dow, tmp$DOW[1])
  
}

axe$PMX <- a.pm
axe$Wind <- a.wind
axe$Speed <- a.speed
axe$TEMP <- a.temp
axe$Rain <- a.rain
axe$WEEK <- a.week
axe$DOW <- a.dow

# Creating county-day level counts

locs <- unique(axe$FIPS)
days <- unique(axe$Date)
alocs <- c()
adays <- c()
axes <- c()
apmx <- c()
awind <- c()
aspeed <- c()
atemp <- c()
arain <- c()
aweek <- c()
adow <- c()

for (l in locs) {
  
  print(l)
  
  tmp <- axe %>% filter(FIPS == l)
  
  for (d in days) {
    
    tmp2 <- tmp %>% filter(Date == d)
    
    alocs <- c(alocs, l)
    adays <- c(adays, d)
    axes <- c(axes, nrow(tmp2))
    apmx <- c(apmx, mean(tmp2$PMX, na.rm = TRUE))
    awind <- c(awind, mean(tmp2$Wind, na.rm = TRUE))
    aspeed <- c(aspeed, mean(tmp2$Speed, na.rm = TRUE))
    atemp <- c(atemp, mean(tmp2$TEMP, na.rm = TRUE))
    arain <- c(arain, mean(tmp2$Rain, na.rm = TRUE))
    aweek <- c(aweek, mean(tmp2$WEEK, na.rm = TRUE))
    adow <- c(adow, mean(tmp2$DOW, na.rm = TRUE))
    
  }
  
}

ax.reg <- as.data.frame(cbind(axes, apmx, awind, aspeed, atemp, arain, aweek, adow))

colnames(ax.reg) <- c('Accidents', 'PMX', 'Wind', 'Speed', 'Temp', 'Rain', 'Week', 'DOW')

ax.reg$FIPS <- alocs
ax.reg$Date <- adays

ax.reg$Icy <- as.integer(ax.reg$Temp < 32) * ax.reg$Rain
ax.reg$Wet <- ax.reg$Rain + 2*ax.reg$Icy
ax.reg <- ax.reg[which(complete.cases(ax.reg)),]

# I should add stops data, too

bab.ybottle.stops <- c()

for (i in 1:nrow(ax.reg)) {
  
  print(i)
  tmp <- data %>% filter(STOP_DATE == ax.reg$Date[i]) %>% filter(FIPS.x == ax.reg$FIPS[i])
  bab.ybottle.stops <- c(bab.ybottle.stops, tmp$Stops[1])

}

ax.reg$Stops <- bab.ybottle.stops

# Creating a month FE

mos <- c()

for (i in 1:nrow(ax.reg)) {
  
  print(i)
  mos <- c(mos, strsplit(ax.reg$Date[i], '/')[[1]][1])
  
}

ax.reg$Month <- mos

# Remove outliers with more than 100 stops/day

ax.reg2 <- ax.reg[which(ax.reg$Stops <= 100),]

# Driver behavior regressions

a1 <- glm(Accidents ~ PMX + log(Stops) + Temp + I(Temp^2) + factor(Wet) + factor(Month) + factor(FIPS) + factor(DOW), data = ax.reg2, family = 'poisson')
a2 <- glm.nb(Accidents ~ PMX + log(Stops) + Temp + I(Temp^2) + factor(Wet) + factor(Month) + factor(FIPS) + factor(DOW), data = ax.reg2)
a3 <- lm(Accidents ~ PMX + log(Stops) + Temp + I(Temp^2) + factor(Wet) + factor(Month) + factor(FIPS) + factor(DOW), data = ax.reg2)
a4 <- ivreg(Accidents ~ PMX + log(Stops) + Temp + I(Temp^2) + factor(Wet) + factor(Month) + factor(FIPS) + factor(DOW) | . - PMX + factor(Wind)*factor(FIPS) + Speed + factor(Week), data = ax.reg2)

a1x <- coeftest(a1, vcov = vcovCL, cluster = ~FIPS)
a2x <- coeftest(a2, vcov = vcovCL, cluster = ~FIPS)
a3x <- coeftest(a3, vcov = vcovCL, cluster = ~FIPS)
a4x <- coeftest(a4, vcov = vcovCL, cluster = ~FIPS)

s1 <- glm(Stops ~ PMX +  Temp + I(Temp^2) + factor(Wet) + factor(Month) + factor(FIPS) + factor(DOW), data = ax.reg2, family = 'poisson')
s2 <- glm.nb(Stops ~ PMX +  Temp + I(Temp^2) + factor(Wet) + factor(Month) + factor(FIPS) + factor(DOW), data = ax.reg2)
s3 <- lm(Stops ~ PMX + Temp + I(Temp^2) + factor(Wet) + factor(Month) + factor(FIPS) + factor(DOW), data = ax.reg2)
s4 <- ivreg(Stops ~ PMX + Temp + I(Temp^2) + factor(Wet) + factor(Month) + factor(FIPS) + factor(DOW) | . - PMX + factor(Wind)*factor(FIPS) + Speed + factor(Week), data = ax.reg2)

s1x <- coeftest(s1, vcov = vcovCL, cluster = ~FIPS)
s2x <- coeftest(s2, vcov = vcovCL, cluster = ~FIPS)
s3x <- coeftest(s3, vcov = vcovCL, cluster = ~FIPS)
s4x <- coeftest(s4, vcov = vcovCL, cluster = ~FIPS)

write.csv(stargazer(a1, a2, a3, a4, s1, s2, s3, s4, type = 'text', omit = c('Month', 'DOW', 'FIPS')), paste0(direc, 'results/robust_unadj.txt'), row.names = FALSE)

write.csv(stargazer(a1x, a2x, a3x, a4x, s1x, s2x, s3x, s4x, type = 'text', omit = c('Month', 'DOW', 'FIPS')), paste0(direc, 'results/robust.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_a <- summary(a4, diagnostics = TRUE)
sum_s <- summary(s4, diagnostics = TRUE)

f_stats <- as.data.frame(c(sum_a[[12]][7], sum_s[[12]][7]))

colnames(f_stats) <- c('F_STAT')

write.csv(f_stats, paste0(direc, 'results/f_stats_robust.txt'), row.names = FALSE)

# Driver behavior regressions with Month x FIPS fixed effects

a1 <- glm(Accidents ~ PMX + log(Stops) + Temp + I(Temp^2) + factor(Wet) + factor(Month)*factor(FIPS) + factor(DOW), data = ax.reg2, family = 'poisson')
a2 <- glm.nb(Accidents ~ PMX + log(Stops) + Temp + I(Temp^2) + factor(Wet) + factor(Month)*factor(FIPS) + factor(DOW), data = ax.reg2)
a3 <- lm(Accidents ~ PMX + log(Stops) + Temp + I(Temp^2) + factor(Wet) + factor(Month)*factor(FIPS) + factor(DOW), data = ax.reg2)
a4 <- ivreg(Accidents ~ PMX + log(Stops) + Temp + I(Temp^2) + factor(Wet) + factor(Month)*factor(FIPS) + factor(DOW) | . - PMX + factor(Wind)*factor(FIPS) + Speed + factor(Week), data = ax.reg2)

a1x <- coeftest(a1, vcov = vcovCL, cluster = ~FIPS)
a2x <- coeftest(a2, vcov = vcovCL, cluster = ~FIPS)
a3x <- coeftest(a3, vcov = vcovCL, cluster = ~FIPS)
a4x <- coeftest(a4, vcov = vcovCL, cluster = ~FIPS)

s1 <- glm(Stops ~ PMX +  Temp + I(Temp^2) + factor(Wet) + factor(Month)*factor(FIPS) + factor(DOW), data = ax.reg2, family = 'poisson')
s2 <- glm.nb(Stops ~ PMX +  Temp + I(Temp^2) + factor(Wet) + factor(Month)*factor(FIPS) + factor(DOW), data = ax.reg2)
s3 <- lm(Stops ~ PMX + Temp + I(Temp^2) + factor(Wet) + factor(Month)*factor(FIPS) + factor(DOW), data = ax.reg2)
s4 <- ivreg(Stops ~ PMX + Temp + I(Temp^2) + factor(Wet) + factor(Month)*factor(FIPS) + factor(DOW) | . - PMX + factor(Wind)*factor(FIPS) + Speed + factor(Week), data = ax.reg2)

s1x <- coeftest(s1, vcov = vcovCL, cluster = ~FIPS)
s2x <- coeftest(s2, vcov = vcovCL, cluster = ~FIPS)
s3x <- coeftest(s3, vcov = vcovCL, cluster = ~FIPS)
s4x <- coeftest(s4, vcov = vcovCL, cluster = ~FIPS)

write.csv(stargazer(a1, a2, a3, a4, s1, s2, s3, s4, type = 'text', omit = c('Month', 'DOW', 'FIPS')), paste0(direc, 'results/robust_unadj_2.txt'), row.names = FALSE)

write.csv(stargazer(a1x, a2x, a3x, a4x, s1x, s2x, s3x, s4x, type = 'text', omit = c('Month', 'DOW', 'FIPS')), paste0(direc, 'results/robust_2.txt'), row.names = FALSE)

# Model diagnostics for IV regressions to get first stage F-statistics

sum_a <- summary(a4, diagnostics = TRUE)
sum_s <- summary(s4, diagnostics = TRUE)

f_stats <- as.data.frame(c(sum_a[[12]][7], sum_s[[12]][7]))

colnames(f_stats) <- c('F_STAT')

write.csv(f_stats, paste0(direc, 'results/f_stats_robust_2.txt'), row.names = FALSE)

# Placebo testing

# Creating data for placebo tests

data2 <- data[,c('TICKETED', 'ARRESTED', 'SEARCHED', 'VEHICLE_SEARCHED', 'PERSON_SEARCHED',
                 'FORCE', 'PMX', 'Stops', 'REASON.FOR.STOP', 'MINOR', 'TEMP', 'Rain', 
                 'RACE', 'ETHNICITY', 'GENDER', 'ENGLISH.SPEAKING', 'RESIDENCY',
                 'AGENCY.NAME', 'JURISDICTION', 'Month', 'WEEK', 'DOW', 'Wind', 'Speed')]

data2 <- data2[which(complete.cases(data2)),]

set.seed(42024)

new.ids <- sample.int(nrow(data2), nrow(data2))

data2$TICKETED_P <- data2$TICKETED[new.ids]
data2$ARRESTED_P <- data2$ARRESTED[new.ids]
data2$SEARCHED_P <- data2$SEARCHED[new.ids]
data2$VEHICLE_SEARCHED_P <- data2$VEHICLE_SEARCHED[new.ids]
data2$PERSON_SEARCHED_P <- data2$PERSON_SEARCHED[new.ids]
data2$FORCE_P <- data2$FORCE[new.ids]

set.seed(3704558)

new.ids2 <- sample.int(nrow(data2), nrow(data2))

data2$PMX_P <- data2$PMX[new.ids2]

# Placebo tests with randomized outcomes

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED_P ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                       + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                       | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED_P ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                         + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                         | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED_P ~ PMX + log(Stops) + log(SEARCHED_P) + factor(REASON.FOR.STOP) + MINOR
                          + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                          + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                          | . - PMX + factor(Wind) + Speed + factor(WEEK) + factor(DOW), data = data2[which(data2$SEARCHED_P == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED_P ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                         + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                         | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED_P ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                        + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                        + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                        | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED_P ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                       + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                       | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE_P ~ PMX + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                      + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                      + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                      | . - PMX + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE_P ~ PMX + log(Stops) + log(ARRESTED_P) + factor(REASON.FOR.STOP) + MINOR
                       + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                       | . - PMX + factor(Wind) + Speed + factor(WEEK) + factor(DOW), data = data2[which(data2$ARRESTED_P == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'Month')),
          paste0(direc, 'results/pm_iv_P1.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, omit = c('AGENCY.NAME', 'JURISDICTION', 'Month')),
          paste0(direc, 'results/pm_iv_tex_P1.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'Month')),
          paste0(direc, 'results/pm_iv_unadj_P1.txt'), row.names = FALSE)

# Placebo tests with randomized pollution

# Running regressions for ticketed as an outcome conditional on being stopped

ticket_iv.xxx <- ivreg(TICKETED ~ PMX_P + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                       + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                       | . - PMX_P + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xticket_iv.xxx <- coeftest(ticket_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being stopped

arrested_iv.xxx <- ivreg(ARRESTED ~ PMX_P + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                         + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                         | . - PMX_P + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xarrested_iv.xxx <- coeftest(arrested_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for arrested as an outcome conditional on being searched

arrestedx_iv.xxx <- ivreg(ARRESTED ~ PMX_P + log(Stops) + log(SEARCHED_P) + factor(REASON.FOR.STOP) + MINOR
                          + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                          + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                          | . - PMX_P + factor(Wind) + Speed + factor(WEEK) + factor(DOW), data = data2[which(data2$SEARCHED_P == 1),])

xarrestedx_iv.xxx <- coeftest(arrestedx_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for being searched in any capacity as an outcome conditional on being stopped

searched_iv.xxx <- ivreg(SEARCHED ~ PMX_P + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                         + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                         + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                         | . - PMX_P + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xsearched_iv.xxx <- coeftest(searched_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the vehicle being searched as an outcome conditional on being stopped

vehicle_iv.xxx <- ivreg(VEHICLE_SEARCHED ~ PMX_P + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                        + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                        + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                        | . - PMX_P + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xvehicle_iv.xxx <- coeftest(vehicle_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the driver being searched as an outcome conditional on being stopped

driver_iv.xxx <- ivreg(PERSON_SEARCHED ~ PMX_P + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                       + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                       | . - PMX_P + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xdriver_iv.xxx <- coeftest(driver_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being stopped

force_iv.xxx <- ivreg(FORCE ~ PMX_P + log(Stops) + factor(REASON.FOR.STOP) + MINOR
                      + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                      + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                      | . - PMX_P + factor(Wind)*factor(JURISDICTION) + Speed + factor(WEEK) + factor(DOW), data = data2)

xforce_iv.xxx <- coeftest(force_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Running regressions for the cop using force as an outcome conditional on being arrested

forcex_iv.xxx <- ivreg(FORCE ~ PMX_P + log(Stops) + log(ARRESTED_P) + factor(REASON.FOR.STOP) + MINOR
                       + TEMP + I(TEMP^2) + Rain + factor(RACE) + factor(ETHNICITY) + factor(GENDER)
                       + factor(ENGLISH.SPEAKING) + factor(RESIDENCY) + factor(AGENCY.NAME)*factor(Month) + factor(JURISDICTION)
                       | . - PMX_P + factor(Wind) + Speed + factor(WEEK) + factor(DOW), data = data2[which(data2$ARRESTED_P == 1),])

xforcex_iv.xxx <- coeftest(forcex_iv.xxx, vcov = vcovCL, cluster = ~AGENCY.NAME)

# Saving results

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'Month')),
          paste0(direc, 'results/pm_iv_P2.txt'), row.names = FALSE)

write.csv(stargazer(xticket_iv.xxx, xarrested_iv.xxx, xarrestedx_iv.xxx, xsearched_iv.xxx, xvehicle_iv.xxx, xdriver_iv.xxx,
                    xforce_iv.xxx, xforcex_iv.xxx, omit = c('AGENCY.NAME', 'JURISDICTION', 'Month')),
          paste0(direc, 'results/pm_iv_tex_P2.txt'), row.names = FALSE)

write.csv(stargazer(ticket_iv.xxx, arrested_iv.xxx, arrestedx_iv.xxx, searched_iv.xxx, vehicle_iv.xxx, driver_iv.xxx,
                    force_iv.xxx, forcex_iv.xxx, type = 'text', omit = c('AGENCY.NAME', 'JURISDICTION', 'Month')),
          paste0(direc, 'results/pm_iv_unadj_P2.txt'), row.names = FALSE)

# Figures for the paper

# Histogram - stops

ggplot(ax.reg2, aes(x = Stops)) +
  theme_bw() +
  theme(axis.line = element_line(color = 'black'),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  ggtitle('Histogram for Stops per County-Day') +
  xlab('Stops') +
  ylab('Frequency') +
  geom_histogram(color = 'black', fill = 'blue4', bins = max(ax.reg2$Stops)+1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(c(0, 1500))

# Histogram - accidents

ggplot(ax.reg2, aes(x = Accidents)) +
  theme_bw() +
  theme(axis.line = element_line(color = 'black'),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  ggtitle('Histogram for Accidents per County-Day') +
  xlab('Accidents') +
  ylab('Frequency') +
  geom_histogram(color = 'black', fill = 'purple4', bins = max(ax.reg2$Accidents)+1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(c(0, 3000))

# Leaflet - the final sample

va$GEOID <- as.integer(va$GEOID)

va$SAMPLE <- ifelse(va$GEOID %in% ax.reg2$FIPS, 1, 0)

p1 <- colorNumeric(palette = c('white', 'red4'), domain = va$SAMPLE)

leaflet(va$geometry) %>%
  addTiles() %>%
  addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = va$SAMPLE, color = 'black', fillColor = p1(va$SAMPLE))

# Leaflet - high PM days per county

hpm <- c()

for (c in unique(va$GEOID)) {
  
  print(c)
  tmp <- ax.reg2 %>% filter(FIPS == c)
  hpm <- c(hpm, sum(tmp$PMX))
  
}

va$HPM <- hpm

p2 <- colorNumeric(palette = 'Reds', domain = va$HPM)

leaflet(va$geometry) %>%
  addTiles() %>%
  addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = va$SAMPLE, color = 'black', fillColor = p2(va$HPM)) %>%
  addLegend(position = 'topleft', pal = p2, values = va$HPM, title = 'High PM Days', opacity = 1)

# Leaflet - agencies per county

ags <- c()

for (c in unique(va$GEOID)) {
  
  print(c)
  tmp <- data %>% filter(FIPS.x == c)
  ags <- c(ags, length(unique(tmp$AGENCY.NAME)))
  
}

va$AGS <- ags

p3 <- colorNumeric(palette = 'Blues', domain = va$AGS)

leaflet(va$geometry) %>%
  addTiles() %>%
  addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = va$SAMPLE, color = 'black', fillColor = p3(va$AGS)) %>%
  addLegend(position = 'topleft', pal = p3, values = va$AGS, title = 'Agencies', opacity = 1)

# Leaflet - Accidents

dents <- c()

for (c in unique(va$GEOID)) {
  
  print(c)
  tmp <- ax.reg2 %>% filter(FIPS == c)
  ifelse(nrow(tmp > 0), dents <- c(dents, mean(tmp$Accidents, na.rm = TRUE)), dents <- c(dents, 0))
  
}

va$DENTS <- dents

p4 <- colorNumeric(palette = 'Oranges', domain = va$DENTS)

leaflet(va$geometry) %>%
  addTiles() %>%
  addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = va$SAMPLE, color = 'black', fillColor = p4(va$DENTS)) %>%
  addLegend(position = 'topleft', pal = p4, values = va$DENTS, title = 'Mean Daily Accidents', opacity = 1)

# Leaflet - Stops

stop <- c()

for (c in unique(va$GEOID)) {
  
  print(c)
  tmp <- ax.reg2 %>% filter(FIPS == c)
  ifelse(nrow(tmp > 0), stop <- c(stop, mean(tmp$Stops, na.rm = TRUE)), stop <- c(stop, 0))
  
}

va$STOP <- stop

p5 <- colorNumeric(palette = 'Purples', domain = va$STOP)

leaflet(va$geometry) %>%
  addTiles() %>%
  addPolygons(weight = 1.0, smoothFactor = 1.0, opacity = 1.0, fillOpacity = va$SAMPLE, color = 'black', fillColor = p5(va$STOP)) %>%
  addLegend(position = 'topleft', pal = p5, values = va$STOP, title = 'Mean Daily Traffic Stops', opacity = 1)

# Summary statistics

data$Asian <- as.integer(data$RACE == 'Asian or Native Hawaiian or Other Pacific Islander')
data$Black <- as.integer(data$RACE == 'Black or African American')
data$Native <- as.integer(data$RACE == 'American Indian or Alaska Native')
data$White <- as.integer(data$RACE == 'White')
data$Hispanic <- as.integer(data$ETHNICITY == 'Hispanic or Latino')
data$Male <- as.integer(data$GENDER == 'Male')

yerp <- c('TICKETED', 'ARRESTED', 'SEARCHED', 'VEHICLE_SEARCHED', 'PERSON_SEARCHED', 'PMX', 'PM',
          'Asian', 'Black', 'Hispanic', 'Native', 'White', 'Male', 'TEMP', 'Rain', 'MINOR')

st <- data[,which(colnames(data) %in% yerp)]
st <- st[complete.cases(st),]
st <- st[,c(6,5,4,3,2,10,8,1,9,11,12,13,14,15,16,7)]

cn <- c('Ticketed', 'Arrested', 'Searched (Any)', 'Driver Searched', 'Vehicle Searched', 'High PM Day', 'PM (ug/m3)',
        'Temperature (F)', 'Rainy Day', 'Asian', 'Black', 'Native American', 'White', 'Hispanic', 'Male', 'Minor')

colnames(st) <- cn

datasummary_skim(st, fmt = '%.3f')

