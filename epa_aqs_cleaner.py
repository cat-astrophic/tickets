# This script cleans raw EPA AQS data

# Importing required modules

import pandas as pd
import os
import re

# Project directory info

filepath = 'D:/EPA/' # Data directory -- update as needed

# Create a list of all folders to explore

folders = os.listdir(filepath + 'raw_data/')

# Main loop

for folder in folders:
    
    # Create a list of all .txt files to parse
    
    files = os.listdir(filepath + 'raw_data/' + folder)
    
    # Initializing some lists
    
    states = []
    counties = []
    sites = []
    pollutants = []
    values = []
    dates = []
    
    # Data work
        
    for file in files:
        
        print(file) # Tracking progress
        f = open(filepath + 'raw_data/' + folder + '/' + file, 'r') # Open file
        d = f.read() # Read in the data as a string
        f.close()
        
        state_ids = [x.start() for x in re.finditer('state_code', d)] # Get positions of all state_codes
        county_ids = [x.start() for x in re.finditer('county_code', d)] # Get positions of all county_codes
        site_ids = [x.start() for x in re.finditer('site_number', d)] # Get positions of all site_numbers
        pollutant_ids = [x.start() for x in re.finditer('parameter_code', d)] # Get positions of all parameter_codes
        value_ids = [x.start() for x in re.finditer('reported_sample_measurement', d)] # Get positions of all reported_measurements
        date_ids = [x.start() for x in re.finditer('sample_begin_date', d)] # Get positions of all sample_dates
        
        states_tmp = [d[idx+14:idx+16] for idx in state_ids]
        counties_tmp = [d[idx+15:idx+18] for idx in county_ids]
        sites_tmp = [d[idx+15:idx+19] for idx in site_ids]
        pollutants_tmp = [d[idx+18:idx+23] for idx in pollutant_ids]
        dates_tmp = [d[idx+21:idx+29] for idx in date_ids]
        values_tmp = []
        
        for idx in value_ids:
            
            a = d[idx+31:].find('"') # The number of digits presented is inconsistent
            values_tmp.append(d[idx+31:idx+31+a])
            
        states = states + states_tmp
        counties = counties + counties_tmp
        sites = sites + sites_tmp
        pollutants = pollutants + pollutants_tmp
        values = values + values_tmp
        dates = dates + dates_tmp
    
    # Create a dataframe containing the complete data set
        
    dates = pd.Series(dates, name = 'Date')
    states = pd.Series(states, name = 'State')
    counties = pd.Series(counties, name = 'County')
    sites = pd.Series(sites, name = 'Site')
    pollutants = pd.Series(pollutants, name = 'Pollutant')
    values = pd.Series(values, name = 'Value')
    df = pd.concat([dates, states, counties, sites, pollutants, values], axis = 1)
    
    # Some of the entries for Value were 'ull,\n      ' whenever the result was null, so we drop those
    
    df = df[df.Value != 'ull,\n      '].reset_index(drop = True)
    
    # Write the dataframe to file
    
    df.to_csv(filepath + 'epa_aqs_data_' + folder + '.csv', index = False)

