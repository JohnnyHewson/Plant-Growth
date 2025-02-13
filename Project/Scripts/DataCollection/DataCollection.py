import csv
import os
import requests
from requests import ConnectTimeout
import pandas as pd
import datetime
import time

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\\..'))

configdir = os.path.join(project_path, 'config.txt')
with open(configdir, 'r') as config_file:
    stages = []
    for line in config_file:
        if line.__contains__(','):
            if line[:line.index(',')] == 'Latitude':
                latitude = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'Longitude':
                longitude = float(line.strip().split(',')[1])

today = datetime.datetime.today() - datetime.timedelta(1)
end_yyyymmdd = int(f"{today.year}{today.month if len(str(today.month)) == 2 else f'0{today.month}'}{today.day if len(str(today.day)) == 2 else f'0{today.day}'}")

###Collect Data###
#Hourly Data (Photosynthetically Active Radiation and Temperature at 2 Metres)
hourly_start_year = 2001
start_time = time.perf_counter()
while True:
    attempt_start = time.perf_counter()
    print(f"Collecting hourly PAR and temperature data at ({latitude},{longitude}) from {hourly_start_year} to now")
    start_date = datetime.datetime(hourly_start_year,1,1)
    start_yyyymmdd = int(f"{start_date.year}{start_date.month if len(str(start_date.month)) == 2 else f'0{start_date.month}'}{start_date.day if len(str(start_date.day)) == 2 else f'0{start_date.day}'}")
    try:
        base_url = r"https://power.larc.nasa.gov/api/temporal/hourly/point?parameters=ALLSKY_SFC_PAR_TOT,T2M&community=RE&longitude={longitude}&latitude={latitude}&start="+f"{start_yyyymmdd.real}&end={end_yyyymmdd.real}"+r"&format=CSV&header=false&time-standard=utc"

        api_request_url = base_url.format(longitude=longitude, latitude=latitude)

        response = requests.get(url=api_request_url, verify=True, timeout=30.00)
        
        print("Response Code:",response)
        if response.ok is False:
            raise ConnectTimeout
    except ConnectTimeout:
        hourly_start_year += 1
        print(f"Connection timed out after {time.perf_counter()-attempt_start} seconds\nRetrying starting from {hourly_start_year}")
    else:
        print(f"Hourly Data Collected in {time.perf_counter()-start_time} seconds")
        break

content = response.content.decode('utf-8').splitlines()
print("Sorting hourly data...")
hourly = pd.DataFrame((row.split(',') for row in content[1:]),columns=content[0].split(','))
hourly['ALLSKY_SFC_PAR_TOT'] = pd.to_numeric(hourly['ALLSKY_SFC_PAR_TOT'], errors='coerce')
hourly['T2M'] = pd.to_numeric(hourly['T2M'], errors='coerce')
hourly['Date'] = hourly['YEAR']+'-'+hourly['MO']+'-'+hourly['DY']
hourly['Date'] = pd.to_datetime(hourly['Date'])
hourly.drop(columns=['YEAR','MO','DY'],inplace=True)
for index,row in hourly.iterrows():
    hourly.loc[index,'JDay'] = row['Date'].timetuple().tm_yday
colnames=list(hourly.columns[:-2].values)
colnames.insert(0,'Date')
colnames.insert(1,'JDay')
hourly = hourly.reindex(columns=colnames)
hourly = hourly.rename(columns={'ALLSKY_SFC_PAR_TOT':'PAR','T2M':'Temp'})

#Daily Data (Daily Minimum Temperautre and Daily Maximum Temperature)
daily_start_year = 1981
start_time = time.perf_counter()
while True:
    attempt_start = time.perf_counter()
    print(f"Collecting daily minimum and maximum temperature data at ({latitude},{longitude}) from {daily_start_year} to now")
    start_date = datetime.datetime(daily_start_year,1,1)
    start_yyyymmdd = int(f"{start_date.year}{start_date.month if len(str(start_date.month)) == 2 else f'0{start_date.month}'}{start_date.day if len(str(start_date.day)) == 2 else f'0{start_date.day}'}")
    try:
        base_url = r"https://power.larc.nasa.gov/api/temporal/daily/point?parameters=T2M_MIN,T2M_MAX&community=RE&longitude={longitude}&latitude={latitude}&start="+f"{start_yyyymmdd.real}&end={end_yyyymmdd.real}"+r"&format=CSV&header=false&time-standard=utc"

        api_request_url = base_url.format(longitude=longitude, latitude=latitude)

        response = requests.get(url=api_request_url, verify=True, timeout=30.00)

        print("Response Code:",response)
        if response.ok is False:
            raise ConnectTimeout
    except ConnectTimeout:
        daily_start_year += 1
        print(f"Connection timed out after {time.perf_counter()-attempt_start} seconds\nRetrying starting from {daily_start_year}")
    else:
        print(f"Daily Data Collected in {time.perf_counter()-start_time} seconds")
        break

content = response.content.decode('utf-8').splitlines()
print("Sorting daily data...")
daily = pd.DataFrame((row.split(',') for row in content[1:]),columns=content[0].split(','))
daily['T2M_MIN'] = pd.to_numeric(daily['T2M_MIN'], errors='coerce')
daily['T2M_MAX'] = pd.to_numeric(daily['T2M_MAX'], errors='coerce')
daily['Date'] = daily['YEAR']+'-'+daily['MO']+'-'+daily['DY']
daily['Date'] = pd.to_datetime(daily['Date'])
daily.drop(columns=['YEAR','MO','DY'],inplace=True)
for index,row in daily.iterrows():
    daily.loc[index,'JDay'] = row['Date'].timetuple().tm_yday
colnames=list(daily.columns[:-2].values)
colnames.insert(0,'Date')
colnames.insert(1,'JDay')
daily = daily.reindex(columns=colnames)
daily = daily.rename(columns={'T2M_MIN':'Min Temp','T2M_MAX':'Max Temp'})

#Test if there is any inconsistencies
i=0
for index,row in daily.iterrows():
    if row['Date'] != hourly.loc[24*i,'Date']:
        continue
    if row['Min Temp'] == min(hourly.loc[24*(i):24*(i+1),'Temp']) or row['Max Temp'] == max(hourly.loc[24*(i):24*(i+1),'Temp']):
        continue
    print(row['Date'])
    print(f"Daily Min Temp: {row['Min Temp']}|||Daily Max Temp: {row['Max Temp']}")
    print(f"Minimum of daily hourly temp: {min(hourly.loc[24*(i):24*(i+1),'Temp'])}|||Max of daily hourly temp: {max(hourly.loc[24*(i):24*(i+1),'Temp'])}")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    i+=1

###Process Data###
#PAR
group = hourly[hourly['PAR'] != float(-999)].groupby(['JDay','HR'], as_index=False)['PAR'].mean()
Average_Hourly_PAR = group.pivot(index='JDay', columns='HR', values='PAR')
Average_Hourly_PAR.columns = Average_Hourly_PAR.columns.astype(int)
Average_Hourly_PAR = Average_Hourly_PAR.sort_index(axis=1)
#Temperature
group = hourly[hourly['Temp'] != float(-999)].groupby(['JDay','HR'], as_index=False)['Temp'].mean()
Average_Hourly_Temp = group.pivot(index='JDay', columns='HR', values='Temp')
Average_Hourly_Temp.columns = Average_Hourly_Temp.columns.astype(int)
Average_Hourly_Temp = Average_Hourly_Temp.sort_index(axis=1)
#Minimum Temperature
Average_Daily_Min_Temp = daily[daily['Min Temp'] != float(-999)].groupby('JDay', as_index=True)['Min Temp'].mean()

#Maximum Temperature
Average_Daily_Max_Temp = daily[daily['Max Temp'] != float(-999)].groupby('JDay', as_index=True)['Max Temp'].mean()

# ###Save Data
# processed_data_path = os.path.join(project_path,'Data','Processed','Average Conditions')
# Average_Hourly_PAR.to_csv(os.path.join(processed_data_path,'Average Hourly PAR'+'.csv'),index=True,header=True)
# Average_Hourly_Temp.to_csv(os.path.join(processed_data_path,'Average Hourly Temperature'+'.csv'),index=True,header=True)
# Average_Daily_Min_Temp.to_csv(os.path.join(processed_data_path,'Average Daily Min Temperature'+'.csv'),index=True,header=True)
# Average_Daily_Max_Temp.to_csv(os.path.join(processed_data_path,'Average Daily Max Temperature'+'.csv'),index=True,header=True)

# print(f"Data saved to {processed_data_path}")