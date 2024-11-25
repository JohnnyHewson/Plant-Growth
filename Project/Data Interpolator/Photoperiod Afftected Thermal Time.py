import pandas as pd
from math import cos, pi, sin, asin, acos, tan
from datetime import datetime, timedelta

# Constants
P_opt = 20  # Optimal photoperiod
P_Base_dict = {'Emergence_to_Double_Ridge': 0, 'Double_Ridge_to_Anthesis': 7}

# Load stage timings and thermal temperature data
stage_timings = pd.read_csv('Stage Timings.csv')
thermal_data = pd.read_csv('ThermalTemperature_78-81_Calculated.csv')

# Convert dates in thermal data to datetime format
thermal_data['Date'] = pd.to_datetime(thermal_data['Date'], format='%Y-%m-%d')

# Ask for the latitude input
Lat = 51.8 #Lat of Rothamsted

def calculate_photoperiod(date):
    """Calculate photoperiod-effective hours (P_H) for a given date."""
    # Calculate Julian day number
    Jday = date.timetuple().tm_yday

    # Declination angle calculation
    theta1 = 2 * pi * (Jday - 80) / 365
    theta2 = 0.0335 * (sin(2 * pi * Jday) - sin(2 * pi * 80))
    theta = theta1 + theta2
    Dec = asin(0.3978 * sin(theta))

    # Photoperiod factor
    D = -0.10453 / (cos(Lat) * cos(Dec))
    P_R = acos(D - tan(Lat) * tan(Dec))
    P_H = 24 * P_R / pi
    return P_H

def calculate_fp(P_H, stage):
    """Calculate FP based on photoperiod hours and growth stage."""
    P_Base = P_Base_dict[stage]
    if P_H > P_Base:
        FP = (P_H - P_Base) / (P_opt - P_Base)
        return FP
    else:
        return 0  # No growth if P_H <= P_Base

# Loop through each plant in stage timings
for plant_id in stage_timings.columns[1:]:
    # Get stage timings for the plant
    emergence_date = pd.to_datetime(stage_timings.loc[stage_timings['Stage'] == 'Emergence', plant_id].values[0])
    double_ridge_date = pd.to_datetime(stage_timings.loc[stage_timings['Stage'] == 'Double Ridge', plant_id].values[0])
    anthesis_date = pd.to_datetime(stage_timings.loc[stage_timings['Stage'] == 'Anthesis', plant_id].values[0])

    # Filter thermal data to relevant dates
    plant_data = thermal_data[(thermal_data['Date'] >= emergence_date) & (thermal_data['Date'] <= anthesis_date)].copy()
    plant_data['FP_Thermal_Time'] = 0  # New column for FP-adjusted thermal time

    # Calculate FP-adjusted thermal time for each day
    for i, row in plant_data.iterrows():
        date = row['Date']
        T_t = row['Thermal Temperature']  # Daily thermal time
        
        # Determine growth stage
        if emergence_date <= date < double_ridge_date:
            stage = 'Emergence_to_Double_Ridge'
        elif double_ridge_date <= date <= anthesis_date:
            stage = 'Double_Ridge_to_Anthesis'
        else:
            continue  # Skip if date is outside relevant range
        
        # Calculate P_H and FP
        P_H = calculate_photoperiod(date)
        FP = calculate_fp(P_H, stage)
        
        # Calculate FP-adjusted thermal time
        plant_data.at[i, 'FP_Thermal_Time'] = T_t * FP

    # Save plant-specific output to a new CSV file
    output_file = f'{plant_id}_FP_Adjusted_Thermal_Time.csv'
    plant_data[['Date', 'FP_Thermal_Time']].to_csv(output_file, index=False)
    print(f"Processed {plant_id} and saved to {output_file}")
