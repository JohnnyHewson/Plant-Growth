﻿from random import seed
import pandas as pd
from math import cos, pi, sin, asin, acos, tan
import datetime

# Constants
P_opt = 20  # Optimal photoperiod-effective hours

# Load the stage timings and thermal temperature CSVs
stage_timings = pd.read_csv('Stage Timings.csv')
thermal_data = pd.read_csv('ThermalTemperature_78-81_Calculated.csv')
min_temp_data = pd.read_csv('Rothamsted Daily Min Temps 01.09.1978-31.08.1981.csv')
max_temp_data = pd.read_csv('Rothamsted Daily Max Temps 01.09.1978-31.08.1981.csv')
thermal_data['Date'] = pd.to_datetime(thermal_data['Date'])  # Ensure dates are in datetime format
min_temp_data['Date'] = pd.to_datetime(min_temp_data['Date'], dayfirst = True)
max_temp_data['Date'] = pd.to_datetime(max_temp_data['Date'], dayfirst = True)

thermal_data = thermal_data.merge(min_temp_data, on='Date', how='outer').merge(max_temp_data, on='Date', how='outer')


# Function to calculate photoperiod-effective hours (P_H)
def calculate_photoperiod(lat, date):
    Jday = date.timetuple().tm_yday
    theta1 = 2 * pi * (Jday - 80) / 365
    theta2 = 0.0335 * (sin(2 * pi * Jday) - sin(2 * pi * 80))
    theta = theta1 + theta2
    Dec = asin(0.3978 * sin(theta))
    D = -0.10453 / (cos(lat * pi / 180) * cos(Dec))
    P_R = acos(D - tan(lat * pi / 180) * tan(Dec))
    return 24 * P_R / pi

# Latitude input (can be modified as needed)
latitude = 51.81 #Rothamsted

# Loop through each plant in the stage timings data
for index, row in stage_timings.iterrows():
    plant_id = row['Plant']
    seeding_date = pd.to_datetime(row['Seeding Date'], format='%d/%m/%Y')
    
    # Calculate the dates for each stage
    emergence_date = seeding_date + pd.Timedelta(days=row['Days till Emergence'])
    double_ridge_date = seeding_date + pd.Timedelta(days=row['Days till Double Ridge'])
    anthesis_date = seeding_date + pd.Timedelta(days=row['Days till Anthesis'])
    maturity_date = seeding_date + pd.Timedelta(days=row['Days till Maturity'])
    
    # Initialize list to store results for this plant
    plant_results = []
    
    # Loop through each day in the thermal data range
    for _, day_data in thermal_data.iterrows():
        date = day_data['Date']
        thermal_temp = day_data['Thermal Temperature']

        ### FP ###
        # Calculate photoperiod-effective hours for the day
        P_H = calculate_photoperiod(latitude, date)
        
        # Determine the FP factor based on the growth stage
        if emergence_date <= date < double_ridge_date:
            P_Base = 0
        elif double_ridge_date <= date < anthesis_date:
            P_Base = 7
        else:
            P_Base = P_opt  # No contribution post-anthesis
        
        # Calculate FP and the photoperiod-affected thermal time
        if P_H >= P_opt:
            FP = 1
        elif emergence_date <= date < anthesis_date:
            FP = min(max((P_H - P_Base) / (P_opt - P_Base), 0), 1)
        else:
            FP = 1
        photoperiod_affected_thermal_time = thermal_temp * FP
        
        ### FV ###
        #Calculate Thermal Hours
        T_min = max(day_data['T_min'],0)
        T_max = max(day_data['T_max'],0)
        V_eff = 0
        V_sat = 33
        V_base = 8

        for r in range(1, 9):
            f_r = (1 / 2) * (1 + cos((90 / 8) * (2 * r - 1) * pi / 180))
            T_H = max(T_min + f_r * (T_max - T_min), 0)
            if -4 <= T_H < 3:
                V_eff += (T_H + 4) / 7
            elif 3 <= T_H <= 10:
                V_eff += 1
            elif 10 < T_H <= 17:
                V_eff += (17 - T_H) / 7
            else:
                V_eff += 0

        if date < seeding_date:
            VDD = 0
        else:
            VDD += (1/8)*V_eff

        if T_max > 30:
            VDD *= 0.5

        if emergence_date <= date < double_ridge_date:
            FV = min(max((VDD-V_base)/(V_sat-V_base), 0), 1)
        else:
            FV = 1

        photovernalized_affected_thermal_time = thermal_temp * FP * FV
        plant_results.append({'Date': date, 'Photovernalized-Affected Thermal Time': photovernalized_affected_thermal_time})

    # Convert results to a DataFrame and save to CSV
    plant_results_df = pd.DataFrame(plant_results)
    output_file = f"Photovernalized Affected Thermal Time/{plant_id}_Photovernalized_Affected_Thermal_Time.csv"
    plant_results_df.to_csv(output_file, index=False)
    print(f"Saved photovernalized-affected thermal time data for {plant_id} to {output_file}")
