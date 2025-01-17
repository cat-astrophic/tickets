# This script analyzes policing data from VA on tickets

# Importing required modules

import pandas as pd
import numpy as np
import requests
import urllib

# Project directory

direc = 'D:/tickets/'

# Reading in the ticket data

data = pd.read_csv(direc + 'data/community-policing-data.csv', error_bad_lines = False)

# Remove observations from out of sample

keeps = [i for i in range(len(data)) if data.STOP_DATE[i][-4:] == '2023']
data = data[data.index.isin(keeps)].reset_index(drop = True)

# Remove observations for drivers under age

data = data[~data.AGE.isin(['Unknown', 'Missing', 'Over 98 Years'])].reset_index(drop = True)
data.AGE = data.AGE.astype(int)
data = data[data.AGE >= 16].reset_index(drop = True)

# Remove observations for passengers

data = data[data['PERSON TYPE'] == 'Driver'].reset_index(drop = True)

# Removing non traffic or equipment violations

data = data[data['REASON FOR STOP'].isin(['Traffic Violation', 'Equipment Violation'])].reset_index(drop = True)

# Adding FIPS to data

counties = ['Accomack CO (001)', 'Albemarle CO (002)', 'Alleghany CO (003)', 'Amelia CO (004)', 'Amherst CO (005)',
            'Appomattox CO (006)',  'Arlington CO (000)', 'Augusta CO (007)', 'Bath CO (008)', 'Bedford CO (009)',
            'Bland CO (010)', 'Botetourt CO (011)', 'Brunswick CO (012)', 'Buchanan CO (013)', 'Buckingham CO (014)',
            'Campbell CO (015)', 'Caroline CO (016)', 'Carroll CO (017)', 'Charles City CO (018)', 'Charlotte CO (019)',
            'Chesterfield CO (020)', 'Clarke CO (021)', 'Craig CO (022)', 'Culpeper CO (023)', 'Cumberland CO (024)',
            'Dickenson CO (025)', 'Dinwiddie CO (026)', 'Essex CO (028)', 'Fairfax CO (029)', 'Fauquier CO (030)', 'Floyd CO (031)',
            'Fluvanna CO (032)', 'Franklin County (033)', 'Frederick CO (034)', 'Giles CO (035)', 'Gloucester CO (036)',
            'Goochland CO (037)', 'Grayson CO (038)', 'Greene CO (039)', 'Greensville CO (040)', 'Halifax CO (041)', 'Hanover CO (042)',
            'Henrico CO (043)', 'Henry CO (044)', 'Highland CO (045)', 'Isle of Wight CO (046)', 'James City CO (047)', 'King and Queen CO (049)', 
            'King George CO (048)', 'King William CO (050)', 'Lancaster CO (051)', 'Lee CO (052)', 'Loudoun CO (053)', 'Louisa CO (054)',
            'Lunenburg CO (055)', 'Madison CO (056)', 'Mathews CO (057)', 'Mecklenburg CO (058)', 'Middlesex CO (059)', 'Montgomery CO (060)',
            'Nelson CO (062)', 'New Kent CO (063)', 'Northampton CO (065)', 'Northumberland CO (066)', 'Nottoway CO (067)', 'Orange CO (068)',
            'Page CO (069)', 'Patrick CO (070)', 'Pittsylvania CO (071)', 'Powhatan CO (072)', 'Prince Edward CO (073)', 'Prince George CO (074)',
            'Prince William CO (076)', 'Pulaski CO (077)', 'Rappahannock CO (078)', 'Richmond CO (079)', 'Roanoke CO (080)', 'Rockbridge CO (081)',
            'Rockingham CO (082)', 'Russell CO (083)', 'Scott CO (084)', 'Shenandoah CO (085)', 'Smyth CO (086)', 'Southampton CO (087)',
            'Spotsylvania CO (088)', 'Stafford CO (089)', 'Surry CO (090)', 'Sussex CO (091)', 'Tazewell CO (092)', 'Warren CO (093)',
            'Washington CO (095)', 'Westmoreland CO (096)', 'Wise CO (097)', 'Wythe CO (098)', 'York CO (099)', 'Alexandria (100)', 'Bristol (101)',
            'Buena Vista     (102)', 'Charlottesville (103)', 'Chesapeake (126)', 'Colonial Heights (105)', 'Danville (107)', 'Falls Church (108)',
            'Franklin City (133)', 'Fredericksburg (109)', 'Galax (110)', 'Hampton (111)', 'Harrisonburg (112)', 'Hopewell (113)',
            'Lynchburg (114)', 'Manassas (138)', 'Martinsville (115)', 'Newport News (116)', 'Norfolk (117)', 'Norton (130)', 'Petersburg (118)',
            'Poquoson (143)', 'Portsmouth (119)', 'Radford (131)', 'Richmond City (120)', 'Roanoke City (121)', 'Salem (137)', 'Stauton (123)',
            'Suffolk (124)', 'Virginia Beach (125)', 'Waynesboro (127)', 'Williamsburg (128)', 'Winchester (129)', 
            'Covington (106)', 'Manassas Park (139)', 'Fairfax City (132)', 'Emporia (135)', 'Lexington (136)']

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

noaa = pd.read_csv(direc + 'data/NOAA/us_data/NOAA_2023.csv')

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

# Making unified data and weather dates as strings

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

pm =  pd.read_csv(direc + 'data/EPA/epa_aqs_data_PM.csv')
pm10 =  pd.read_csv(direc + 'data/EPA/epa_aqs_data_PM10.csv')
o3 =  pd.read_csv(direc + 'data/EPA/epa_aqs_data_O3.csv')
no2 =  pd.read_csv(direc + 'data/EPA/epa_aqs_data_NO2.csv')
co =  pd.read_csv(direc + 'data/EPA/epa_aqs_data_CO.csv')

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

# Prepping pollution data

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

# Adding pollution to data

data = pd.concat([data, pd.Series([int(d) - 10000 for d in data.Merge_Date], name = 'Merge_Date2')], axis = 1)

data = pd.concat([data, pd.Series([str(data.Merge_Date[i]) + str(data.FIPS[i]) for i in range(len(data))], name = 'AAAAAAAAAAAAAA')], axis = 1)
data = pd.concat([data, pd.Series([str(data.Merge_Date2[i]) + str(data.FIPS[i]) for i in range(len(data))], name = 'BBBBBBBBBBBBBB')], axis = 1)

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

pm_py = []
pm10_py = []
o3_py = []
no2_py = []
co_py = []

for i in range(len(data)):
    
    print(i)
    
    tmp_pm = pm[pm.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]]
    tmp_pm10 = pm10[pm10.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]]
    tmp_o3 = o3[o3.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]]
    tmp_no2 = no2[no2.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]]
    tmp_co = co[co.AAAAAAAAAAAAAA == data.AAAAAAAAAAAAAA[i]]
    
    tmp_pm2 = pm[pm.AAAAAAAAAAAAAA == data.BBBBBBBBBBBBBB[i]]
    tmp_pm102 = pm10[pm10.AAAAAAAAAAAAAA == data.BBBBBBBBBBBBBB[i]]
    tmp_o32 = o3[o3.AAAAAAAAAAAAAA == data.BBBBBBBBBBBBBB[i]]
    tmp_no22 = no2[no2.AAAAAAAAAAAAAA == data.BBBBBBBBBBBBBB[i]]
    tmp_co2 = co[co.AAAAAAAAAAAAAA == data.BBBBBBBBBBBBBB[i]]
    
    if len(tmp_pm) > 0:
        
        pm_xxx.append(np.mean(tmp_pm.PM))
        
    else:
        
        pm_xxx.append(None)
        
    if len(tmp_pm10) > 0:
        
        pm10_xxx.append(np.mean(tmp_pm10.PM10))
        
    else:
        
        pm10_xxx.append(None)
        
    if len(tmp_o3) > 0:
        
        o3_xxx.append(np.mean(tmp_o3.O3))
        
    else:
        
        o3_xxx.append(None)
        
    if len(tmp_no2) > 0:
        
        no2_xxx.append(np.mean(tmp_no2.NO2))
        
    else:
        
        no2_xxx.append(None)
        
    if len(tmp_co) > 0:
        
        co_xxx.append(np.mean(tmp_co.CO))
        
    else:
        
        co_xxx.append(None)
        
    if len(tmp_pm2) > 0:
        
        pm_py.append(np.mean(tmp_pm2.PM))
        
    else:
        
        pm_py.append(None)
        
    if len(tmp_pm102) > 0:
        
        pm10_py.append(np.mean(tmp_pm102.PM10))
        
    else:
        
        pm10_py.append(None)
        
    if len(tmp_o32) > 0:
        
        o3_py.append(np.mean(tmp_o32.O3))
        
    else:
        
        o3_py.append(None)
        
    if len(tmp_no22) > 0:
        
        no2_py.append(np.mean(tmp_no22.NO2))
        
    else:
        
        no2_py.append(None)
        
    if len(tmp_co2) > 0:
        
        co_py.append(np.mean(tmp_co2.CO))
        
    else:
        
        co_py.append(None)

pm_xxx = pd.Series(pm_xxx, name = 'PM')
pm10_xxx = pd.Series(pm10_xxx, name = 'PM10')
o3_xxx = pd.Series(o3_xxx, name = 'Ozone')
no2_xxx = pd.Series(no2_xxx, name = 'NO2')
co_xxx = pd.Series(co_xxx, name = 'CO')

pm_py = pd.Series(pm_py, name = 'PM_PY')
pm10_py = pd.Series(pm10_py, name = 'PM10_PY')
o3_py = pd.Series(o3_py, name = 'Ozone_PY')
no2_py = pd.Series(no2_py, name = 'NO2_PY')
co_py = pd.Series(co_py, name = 'CO_PY')

data = pd.concat([data, pm_xxx, pm10_xxx, o3_xxx, no2_xxx, co_xxx, pm_py, pm10_py, o3_py, no2_py, co_py], axis = 1)

# Adding a month and year variable

month = [x[:x.find('/')] for x in data.STOP_DATE]
year = [x[-4:] for x in data.STOP_DATE]

month = pd.Series(month, name = 'MONTH')
year = pd.Series(year, name = 'YEAR')

data = pd.concat([data, month, year], axis = 1)

# Saving data

data.to_csv(direc + 'data/data.csv', index = False)

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
jtemp = []
jmaxt = []
jmint = []
jprcp = []
jprcpb = []
jpm = []
jpm10 = []
jo3 = []
jco = []
jno2 = []
jpmpy = []
jpm10py = []
jo3py = []
jcopy = []
jno2py = []

for j in tmp.JURISDICTION.unique():
    
    print(j)
    
    jtmp = tmp[tmp.JURISDICTION == j]
    
    ags_x = jtmp['AGENCY NAME'].unique()
    
    for ag in ags_x:
        
        jjtmp = jtmp[jtmp['AGENCY NAME'] == ag].reset_index(drop = True)
        
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
            jtemp.append(jjtmp.TEMP[0])
            jmaxt.append(jjtmp.MAX[0])
            jmint.append(jjtmp.MIN[0])
            jprcp.append(jjtmp.PRCP[0])
            jprcpb.append(int(jjtmp.PRCP[0] > 0))
            jpm.append(jjtmp.PM[0])
            jpm10.append(jjtmp.PM10[0])
            jo3.append(jjtmp.Ozone[0])
            jco.append(jjtmp.CO[0])
            jno2.append(jjtmp.NO2[0])
            jpmpy.append(jjtmp.PM_PY[0])
            jpm10py.append(jjtmp.PM10_PY[0])
            jo3py.append(jjtmp.Ozone_PY[0])
            jcopy.append(jjtmp.CO_PY[0])
            jno2py.append(jjtmp.NO2_PY[0])

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
jtemp = pd.Series(jtemp, name = 'Temperature')
jmaxt = pd.Series(jmaxt, name = 'Max_Temp')
jmint = pd.Series(jmint, name = 'Min_Temp')
jprcp = pd.Series(jprcp, name = 'Precipitation')
jprcpb = pd.Series(jprcpb, name = 'Precipitation_Bin')
jpm = pd.Series(jpm, name = 'PM')
jpm10 = pd.Series(jpm10, name = 'PM10')
jo3 = pd.Series(jo3, name = 'Ozone')
jco = pd.Series(jco, name = 'CO')
jno2 = pd.Series(jno2, name = 'NO2')
jpmpy = pd.Series(jpmpy, name = 'PM_PY')
jpm10py = pd.Series(jpm10py, name = 'PM10_PY')
jo3py = pd.Series(jo3py, name = 'Ozone_PY')
jcopy = pd.Series(jcopy, name = 'CO_PY')
jno2py = pd.Series(jno2py, name = 'NO2_PY')

counts = pd.concat([jur, ags, dat, stops, tix, veh, pers, arr, force, fox,
                    jtemp, jmaxt, jmint, jprcp, jprcpb, jpm, jpm10, jo3,
                    jco, jno2, jpmpy, jpm10py, jo3py, jcopy, jno2py], axis = 1)

# Save the counts data

counts.to_csv(direc + 'data/counts_data.csv', index = False)

