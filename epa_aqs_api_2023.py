# This scripts pulls data from the EPA AQS API

# Importing required modules

import time
import urllib
from bs4 import BeautifulSoup as bs

# Defining some variables

filepath = 'D:/tickets/data/EPA/raw_data/'
email = '' # Username for EPA AQS API
key = '' # User key
params = ['88101', '42602', '44201', '81102', '42101'] # Desired pollutants
param_names = ['PM', 'NO2', 'O3', 'PM10', 'CO'] # Names for these pollutants
months = [str(i) if i>9 else '0'+str(i) for i in range(1,13)] # Months
edates = ['31', '28', '31', '30', '31', '30', '31', '31', '30', '31', '30', '31'] # End dates
states = ['51'] # States FIPS codes

# Defining the url components

url0 = 'https://aqs.epa.gov/data/api/transactionsSample/byState?email='
url1 = '&key='
url2 = '&param='
url3 = '&bdate='
url4 = '&edate='
url5 = '&state='

# Retrieving the data

for param in params:
    
    missed = [] # Stores pollution-state-months which do not come through the API
    
    for state in states:
        
        for year in range(2022,2024):
            
            for month in months:
                
                try: # For certain pollutant-state-months the data does not pull through and must be collected manually
                    
                    print('Pollutant / State / Year / Month :: ' + param_names[params.index(param)] + ' / ' + str(states.index(state)+1) + ' / ' + str(year) + ' / ' + month) # Status update
                    bdate = str(year) + month + '01' # Define bdate
                    edate = str(year) + month + edates[months.index(month)] # Define edate
                    url = url0 + email + url1 + key + url2 + param + url3 + bdate + url4 + edate + url5 + state # Define url
                    page = urllib.request.Request(url, headers = {'User-Agent': 'Mozilla/5.0'}) # Make request
                    response = urllib.request.urlopen(page) # Request data
                    soup = bs(response, 'html.parser') # Clean data
                    
                    with open(filepath + '/' + param_names[params.index(param)] + '/' + param_names[params.index(param)] + '__' + state + '__' + str(year) + '__' + month + '.txt', 'a') as file:
                        
                        file.write(str(soup)) # Write raw data to file
                        file.close() # Close file
                        
                except:
                    
                    missed.append(str(param_names[params.index(param)]) + ' :: ' + str(state) + ' :: ' + str(month) + '/' + str(year))
                    
                    continue
                
                time.sleep(1) # Be kind to your local API

