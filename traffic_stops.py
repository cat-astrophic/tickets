# This script analyzes policing data from VA on tickets

# Importing required modules

import pandas as pd
from datetime import datetime
import requests
import urllib

# Project directory

direc = 'D:/tickets/'

# Reading in the ticket data

data = pd.read_csv(direc + 'data/community-policing-data.csv')

# Remove observations from out of sample

drops_early = [i for i in range(len(data)) if datetime.strptime(data.STOP_DATE[i], '%m/%d/%Y') < datetime(2022,12,31,0,0)]
drops_late = [i for i in range(len(data)) if datetime.strptime(data.STOP_DATE[i], '%m/%d/%Y') >= datetime(2024,1,1,0,0)]
drops = list(set(drops_early + drops_late))
keeps = [i for i in range(len(data)) if i not in drops]
data = data[data.index.isin(keeps)].reset_index(drop = True)

# Remove observations for drivers under age

data = data[data.AGE >= 16].reset_index(drop = True)

# Remove observations for passengers

data = data[data['PERSON TYPE'] == 'DRIVER'].reset_index(drop = True)

# Removing non traffic or equipment violations

data = data[data['REASON FOR STOP'].isin(['TRAFFIC VIOLATION', 'EQUIPMENT VIOLATION'])].reset_index(drop = True)

# Removing observations with non-text jurisdiction

these = list(sorted(data.JURISDICTION.unique()))[155:]

data = data[data.JURISDICTION.isin(these)].reset_index(drop = True)

# Adding FIPS to data

counties = ['ACCOMACK CO', 'ALBEMARLE CO', 'ALLEGHANY CO', 'AMELIA CO', 'AMHERST CO', 'APPOMATTOX CO',
            'ARLINGTON CO', 'AUGUSTA CO', 'BATH CO', 'BEDFORD CO', 'BLAND CO', 'BOTETOURT CO',
            'BRUNSWICK CO', 'BUCHANAN CO', 'BUCKINGHAM CO', 'CAMPBELL CO', 'CAROLINE CO', 'CARROLL CO',
            'CHARLES CITY CO', 'CHARLOTTE CO', 'CHESTERFIELD CO', 'CLARKE CO', 'CRAIG CO', 'CULPEPER CO',
            'CUMBERLAND CO', 'DICKENSON CO', 'DINWIDDIE CO', 'ESSEX CO', 'FAIRFAX CO', 'FAUQUIER CO',
            'FLOYD CO', 'FLUVANNA CO', 'FRANKLIN COUNTY', 'FREDERICK CO', 'GILES CO', 'GLOUCESTER CO',
            'GOOCHLAND CO', 'GRAYSON CO', 'GREENE CO', 'GREENSVILLE CO', 'HALIFAX CO', 'HANOVER CO',
            'HENRICO CO', 'HENRY CO', 'HIGHLAND CO', 'ISLE OF WIGHT CO', 'JAMES CITY CO', 'KING AND QUEEN CO',
            'KING GEORGE CO', 'KING WILLIAM CO', 'LANCASTER CO', 'LEE CO', 'LOUDOUN CO', 'LOUISA CO',
            'LUNENBURG CO', 'MADISON CO', 'MATHEWS CO', 'MECKLENBURG CO', 'MIDDLESEX CO', 'MONTGOMERY CO',
            'NELSON CO', 'NEW KENT CO', 'NORTHAMPTON CO', 'NORTHUMBERLAND CO', 'NOTTOWAY CO', 'ORANGE CO',
            'PAGE CO', 'PATRICK CO', 'PITTSYLVANIA CO', 'POWHATAN CO', 'PRINCE EDWARD CO', 'PRINCE GEORGE CO',
            'PRINCE WILLIAM CO', 'PULASKI CO', 'RAPPAHANNOCK CO', 'RICHMOND CO', 'ROANOKE CO', 'ROCKBRIDGE CO',
            'ROCKINGHAM CO', 'RUSSELL CO', 'SCOTT CO', 'SHENANDOAH CO', 'SMYTH CO', 'SOUTHAMPTON CO',
            'SPOTSYLVANIA CO', 'STAFFORD CO', 'SURRY CO', 'SUSSEX CO', 'TAZEWELL CO', 'WARREN CO',
            'WASHINGTON CO', 'WESTMORELAND CO', 'WISE CO', 'WYTHE CO', 'YORK CO', 'ALEXANDRIA', 'BRISTOL',
            'BUENA VISTA', 'CHARLOTTESVILLE', 'CHESAPEAKE', 'COLONIAL HEIGHTS', 'DANVILLE', 
            'FALLS CHURCH', 'FRANKLIN CITY', 'FREDERICKSBURG', 'GALAX', 'HAMPTON', 'HARRISONBURG',
            'HOPEWELL', 'LYNCHBURG', 'MANASSAS', 'MARTINSVILLE', 'NEWPORT NEWS', 'NORFOLK', 'NORTON',
            'PETERSBURG', 'POQUOSON', 'PORTSMOUTH', 'RADFORD', 'RICHMOND CITY', 'ROANOKE CITY', 'SALEM',
            'STAUNTON', 'SUFFOLK', 'VIRGINIA BEACH', 'WAYNESBORO', 'WILLIAMSBURG', 'WINCHESTER',
            'COVINGTON', 'MANASSAS PARK', 'FAIRFAX CITY', 'EMPORIA', 'LEXINGTON']

fips_codes = ['51001', '51003', '51005', '51007', '51009', '51011', '51013', '51015', '51017', '51019',
              '51021', '51023', '51025', '51027', '51029', '51031', '51033', '51035', '51036', '51037',
              '51041', '51043', '51045', '51047', '51049', '51051', '51053', '51057', '51059', '51061',
              '51063', '51065', '51067', '51069', '51071', '51073', '51075', '51077', '51079', '51081',
              '51083', '51085', '51087', '51089', '51091', '51093', '51095', '51097', '51099', '51101',
              '51103', '51105', '51107', '51109', '51111', '51113', '51115', '51117', '51119', '51121',
              '51125', '51127', '51131', '51133', '51135', '51137', '51139', '51141', '51143', '51145',
              '51147', '51149', '51153', '51155', '51157', '51159', '51161', '51163', '51165', '51167',
              '51169', '51171', '51173', '51175', '51177', '51179', '51181', '51183', '51185', '51187',
              '51191', '51193', '51195', '51197', '51199', '51510', '51520', '51530', '51540', '51550',
              '51570', '51590', '51610', '51620', '51630', '51640', '51650', '51660', '51670', '51680',
              '51683', '51690', '51700', '51710', '51720', '51730', '51735', '51740', '51750', '51760',
              '51770', '51775', '51790', '51800', '51810', '51820', '51830', '51840',
              '51580', '51685', '51600', '51081', '51678']

fips = [fips_codes[counties.index(j)] for j in data.JURISDICTION]

data = pd.concat([data, pd.Series(fips, name = 'FIPS')], axis = 1)

# Bringing in weather data

noaa = pd.read_csv('D:/NOAA/us_data/NOAA_2023.csv')

# Adding FIPS codes to noaa

nerds = [str(noaa.LATITUDE[i]) + str(noaa.LONGITUDE[i]) for i in range(len(noaa))]
ner_df = pd.concat([pd.Series(nerds), pd.Series([i for i in range(len(nerds))])], axis = 1)
ner_df.columns = ['nerds', 'candy']
ner_df = ner_df.drop_duplicates(['nerds']).reset_index(drop = True)

noaa_fips = []

for corn in ner_df.candy:
    
    print(list(ner_df.candy).index(corn))
    params = urllib.parse.urlencode({'latitude': noaa.LATITUDE[corn], 'longitude': noaa.LONGITUDE[corn], 'format': 'json'})
    url = 'https://geo.fcc.gov/api/census/block/find?' + params
    response = requests.get(url)
    content = response.json()
    noaa_fips.append(content['County']['FIPS'])

fips = [noaa_fips[list(ner_df.nerds).index(nerd)] for nerd in nerds]

noaa = pd.concat([noaa, pd.Series(fips, name = 'FIPS')], axis = 1)

# Making unified data and weather dates as strings bc hazy

def data_dates(inp):
    
    y = inp.index('/')
    
    mon = inp[:y]
    day = inp[-1*len(inp) + y + 1:-5]
    
    if len(mon) == 1:
        
        mon = '0' + mon
        
    if len(day) == 1:
        
        day = '0' + day
       
    xxx = inp[-4:] + mon + day
    
    return xxx

data_xxx = [data_dates(d) for d in data.STOP_DATE]
noaa_xxx = [d[:4] + d[5:7] + d[8:] for d in noaa.DATE]

data = pd.concat([data, pd.Series(data_xxx, name = 'Merge_Date')], axis = 1)
noaa = pd.concat([noaa, pd.Series(noaa_xxx, name = 'Merge_Date')], axis = 1)

# Adding weather to data

noaa = noaa.drop_duplicates(['Merge_Date', 'FIPS']).reset_index(drop = True)
data = pd.merge(data, noaa, on = ['Merge_Date', 'FIPS'], how = 'left')

# Bringing in pollution data

pm =  pd.read_csv(direc + 'data/epa_aqs_data_PM.csv')
pm10 =  pd.read_csv(direc + 'data/epa_aqs_data_PM10.csv')
o3 =  pd.read_csv(direc + 'data/epa_aqs_data_O3.csv')
no2 =  pd.read_csv(direc + 'data/epa_aqs_data_NO2.csv')
co =  pd.read_csv(direc + 'data/epa_aqs_data_CO.csv')

# Adding FIPS to pollution data

def fips_fix(inp):
    
    if len(str(inp)) == 1:
        
        fip = '51' + '00' + str(inp)
        
    elif len(str(inp)) == 2:
        
        fip = '51' + '0' + str(inp)
        
    else:
        
        fip = '51' + str(inp)
    
    return fip

pm = pd.concat([pm, pd.Series([fips_fix(pm.County[i]) for i in range(len(pm))], name = 'FIPS')], axis = 1)
pm10 = pd.concat([pm10, pd.Series([fips_fix(pm10.County[i]) for i in range(len(pm10))], name = 'FIPS')], axis = 1)
o3 = pd.concat([o3, pd.Series([fips_fix(o3.County[i]) for i in range(len(o3))], name = 'FIPS')], axis = 1)
no2 = pd.concat([no2, pd.Series([fips_fix(no2.County[i]) for i in range(len(no2))], name = 'FIPS')], axis = 1)
co = pd.concat([co, pd.Series([fips_fix(co.County[i]) for i in range(len(co))], name = 'FIPS')], axis = 1)

# Adding pollution to data

pm['Date'] = pm['Date'].astype(str)
pm10['Date'] = pm10['Date'].astype(str)
o3['Date'] = o3['Date'].astype(str)
no2['Date'] = no2['Date'].astype(str)
co['Date'] = co['Date'].astype(str)

pm = pm[['Date', 'FIPS', 'Value']]
pm10 = pm10[['Date', 'FIPS', 'Value']]
o3 = o3[['Date', 'FIPS', 'Value']]
no2 = no2[['Date', 'FIPS', 'Value']]
co = co[['Date', 'FIPS', 'Value']]

pm.columns = ['Date', 'FIPS', 'PM']
pm10.columns = ['Date', 'FIPS', 'PM10']
o3.columns = ['Date', 'FIPS', 'O3']
no2.columns = ['Date', 'FIPS', 'NO2']
co.columns = ['Date', 'FIPS', 'CO']

pm = pm[pm.PM > 0].reset_index(drop = True)
pm10 = pm10[pm10.PM10 > 0].reset_index(drop = True)
o3 = o3[o3.O3 > 0].reset_index(drop = True)
no2 = no2[no2.NO2 > 0].reset_index(drop = True)
co = co[co.CO > 0].reset_index(drop = True)

data = pd.concat([data, pd.Series([i+1 for i in range(len(data))], name = 'ID')], axis = 1)
data = pd.concat([data, pd.Series([str(data.Merge_Date[i]) + str(data.FIPS[i]) for i in range(len(data))], name = 'AAAAAAAAAAAAAA')], axis = 1)

pm = pd.concat([pm, pd.Series([str(pm.Date[i]) + str(pm.FIPS[i]) for i in range(len(pm))], name = 'AAAAAAAAAAAAAA')], axis = 1)
pm10 = pd.concat([pm10, pd.Series([str(pm10.Date[i]) + str(pm10.FIPS[i]) for i in range(len(pm10))], name = 'AAAAAAAAAAAAAA')], axis = 1)
o3 = pd.concat([o3, pd.Series([str(o3.Date[i]) + str(o3.FIPS[i]) for i in range(len(o3))], name = 'AAAAAAAAAAAAAA')], axis = 1)
no2 = pd.concat([no2, pd.Series([str(no2.Date[i]) + str(no2.FIPS[i]) for i in range(len(no2))], name = 'AAAAAAAAAAAAAA')], axis = 1)
co = pd.concat([co, pd.Series([str(co.Date[i]) + str(co.FIPS[i]) for i in range(len(co))], name = 'AAAAAAAAAAAAAA')], axis = 1)

pm_xxx = []
pm10_xxx = []
o3_xxx = []
no2_xxx = []
co_xxx = []

for i in range(len(data)):
    
    print(i)
    
    tmp_pm = pm[pm.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]].reset_index(drop = True)
    tmp_pm10 = pm10[pm10.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]].reset_index(drop = True)
    tmp_o3 = o3[o3.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]].reset_index(drop = True)
    tmp_no2 = no2[no2.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]].reset_index(drop = True)
    tmp_co = co[co.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]].reset_index(drop = True)
    
    if len(tmp_pm) > 0:
        
        pm_xxx.append(tmp_pm.PM[0])
        
    else:
        
        pm_xxx.append(None)
        
    if len(tmp_pm10) > 0:
        
        pm10_xxx.append(tmp_pm10.PM10[0])
        
    else:
        
        pm10_xxx.append(None)
        
    if len(tmp_o3) > 0:
        
        o3_xxx.append(tmp_o3.O3[0])
        
    else:
        
        o3_xxx.append(None)
        
    if len(tmp_no2) > 0:
        
        no2_xxx.append(tmp_no2.NO2[0])
        
    else:
        
        no2_xxx.append(None)
        
    if len(tmp_co) > 0:
        
        co_xxx.append(tmp_co.CO[0])
        
    else:
        
        co_xxx.append(None)

pm_xxx = pd.Series(pm_xxx, name = 'PM')
pm10_xxx = pd.Series(pm10_xxx, name = 'PM10')
o3_xxx = pd.Series(o3_xxx, name = 'Ozone')
no2_xxx = pd.Series(no2_xxx, name = 'NO2')
co_xxx = pd.Series(co_xxx, name = 'CO')

data = pd.concat([data, pm_xxx, pm10_xxx, o3_xxx, no2_xxx, co_xxx], axis = 1)

# Adding a month and year variable

month = [x[:x.find('/')] for x in data.STOP_DATE]
year = [x[-4:] for x in data.STOP_DATE]

month = pd.Series(month, name = 'MONTH')
year = pd.Series(year, name = 'YEAR')

data = pd.concat([data, month, year], axis = 1)

# Saving data

data.to_csv(direc + 'data/data_raw.csv', index = False)

# Creating daily agency-by-jurisdiction stops data

tmp = data[data.PM > 0].reset_index(drop = True)

days = list(set(tmp.DATE))[1:]

jur = []
ags = []
dat = []
stops = []
tix = []
sear = []
veh = []
pers = []
arr = []
force = []
fox = []

for j in tmp.JURISDICTION.unique():
    
    print(j)
    
    jtmp = tmp[tmp.JURISDICTION == j]
    
    ags_x = jtmp['AGENCY NAME'].unique()
    
    for ag in ags_x:
        
        jjtmp = jtmp[jtmp['AGENCY NAME'] == ag]
        
        for d in days:
            
            jjtmpx = jjtmp[jjtmp.DATE == d]
            tixtmp = jjtmpx[jjtmpx['ACTION TAKEN'] == 'CITATION/SUMMONS']
            vstmp = jjtmpx[jjtmpx['VEHICLE SEARCHED'] == 'YES']
            dstmp = jjtmpx[jjtmpx['PERSON SEARCHED'] == 'YES']
            arrtmp = jjtmpx[jjtmpx['ACTION TAKEN'] == 'ARREST']
            fortmp = jjtmpx[jjtmpx['FORCE USED BY OFFICER'] == 'YES']
            foxtmp = arrtmp[arrtmp['PERSON SEARCHED'] == 'YES']
            
            jur.append(j)
            ags.append(ag)
            dat.append(d)
            stops.append(len(jjtmpx))
            tix.append(len(tixtmp))
            sear.append(len(vstmp) + len(dstmp))
            veh.append(len(vstmp))
            pers.append(len(dstmp))
            arr.append(len(arrtmp))
            force.append(len(fortmp))
            fox.append(len(foxtmp))

jur = pd.Series(jur, name = 'Jurisdiction')
ags = pd.Series(ags, name = 'Agency')
dat = pd.Series(dat, name = 'Date')
stops = pd.Series(stops, name = 'Stops')
tix = pd.Series(tix, name = 'Tickets')
sear = pd.Series(sear, name = 'Searches')
veh = pd.Series(veh, name = 'Vehicles_Searched')
pers = pd.Series(pers, name = 'Persons_Searched')
arr = pd.Series(arr, name = 'Arrests')
force = pd.Series(force, name = 'Force_Used')
fox = pd.Series(fox, name = 'Force_Used_Arrested')

counts = pd.concat([jur, ags, dat, stops, tix, veh, pers, arr, force, fox], axis = 1)

# Save the counts data

counts.to_csv(direc + 'data/counts_data.csv', index = False)

