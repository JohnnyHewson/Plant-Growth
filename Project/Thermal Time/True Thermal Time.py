from random import seed
import pandas as pd
from math import cos, pi, sin, asin, acos, tan
import datetime
import os

# Load the stage timings and temperature data
stage_timings = pd.read_csv('Stage Timings.csv')

min_temp_data = pd.read_csv('Rothamsted Daily Min Temps 01.09.1978-31.08.1981.csv')
max_temp_data = pd.read_csv('Rothamsted Daily Max Temps 01.09.1978-31.08.1981.csv')

min_temp_data['Date'] = pd.to_datetime(min_temp_data['Date'], dayfirst=True)
max_temp_data['Date'] = pd.to_datetime(max_temp_data['Date'], dayfirst=True)

thermal_data = max_temp_data.merge(min_temp_data, on='Date', how='outer')

# Load stage details from config.txt
stages = []
configdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config.txt'))
with open(configdir, 'r') as config_file:
    for line in config_file:
        if line.__contains__(','):
            if line[:line.index(',')] in ['Seeding',
                                          'Emergence',
                                          'Double Ridge',
                                          'Anthesis',
                                          'Maturity']:
                stage_name, T_base, degree_days = line.strip().split(',')
                stages.append({
                    'Stage': stage_name,
                    'T_Base': float(T_base),
                    'Degree_Days': float(degree_days)
                    })


# Constants
P_opt = 20  # Optimal photoperiod-effective hours
V_sat = 33
V_base = 8

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

# Function to calculate daily thermal time (degree days)
def calculate_thermal_time(T_min, T_max, T_base):
    T_opt = 26
    TD_max = 37
    T_t = 0
    for r in range(1, 9):
        f_r = (1 / 2) * (1 + cos((90 / 8) * (2 * r - 1) * pi / 180))
        T_H = max(T_min + f_r * (T_max - T_min), 0)  # Degree Celsius
        if T_H < T_opt:
            T_t += T_H - T_base
        elif T_H == T_opt:
            T_t += T_opt - T_base
        elif T_H < TD_max:
            T_t += (T_opt - T_base) * (TD_max - T_H) / (TD_max - T_opt)
        else:
            T_t += 0
    return max((1 / 8) * T_t, 0)  # Degree Celsius Days

# Latitude input (can be modified as needed)
latitude = 51.81  # Rothamsted

# Loop through each plant in the stage timings data
for index, row in stage_timings.iterrows():
    plant_id = row['Plant']
    seeding_date = pd.to_datetime(row['Seeding Date'], format='%d/%m/%Y')

    # Determine the growing year range
    year_start = datetime.datetime(seeding_date.year, 9, 1)
    year_end = datetime.datetime(seeding_date.year + 1, 8, 31)

    # Initialize cumulative degree days and stage tracking
    cumulative_PVTt = 0
    VDD = 0
    current_stage_index = 0
    current_stage = stages[current_stage_index]

    # Filter thermal data to only include rows within the growing year
    filtered_thermal_data = thermal_data[
    (thermal_data['Date'] >= seeding_date) & (thermal_data['Date'] <= year_end)
    ]

    # Initialize list to store results for this plant
    plant_results = []
    stage_days = 0
    for _, day_data in filtered_thermal_data.iterrows():
        date = day_data['Date']
        T_min = max(day_data['T_min'], 0)
        T_max = max(day_data['T_max'], 0)
        T_base = current_stage['T_Base']
        # Calculate daily thermal time (degree days)
        thermal_time = calculate_thermal_time(T_min, T_max, T_base)

        # Calculate photoperiod (P_H)
        P_H = calculate_photoperiod(latitude, date)

        if current_stage_index == 1:
            P_Base = 0
        elif current_stage_index == 2:
            P_Base = 7
        else:
            P_Base = P_opt  # No contribution post-anthesis
        
        # Calculate FP and the photoperiod-affected thermal time
        if P_H >= P_opt:
            FP = 1
        elif 1 <= current_stage_index < 3:
            FP = min(max((P_H - P_Base) / (P_opt - P_Base), 0), 1)
        else:
            FP = 1

        # Calculate vernalization (VDD)
        V_eff = 0
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
        
        VDD += (1/8)*V_eff

        if T_max > 30:
            VDD *= 0.5

        if current_stage_index > 0:
            FV = min(max((VDD-V_base)/(V_sat-V_base), 0), 1)
        else:
            FV = 1

        PVTt = thermal_time * FP * FV
        cumulative_PVTt += PVTt
        stage_days += 1

        # Store the daily results
        plant_results.append({
            'Date': date,
            'Stage': current_stage['Stage'],
            'Stage Length': stage_days,
            'Daily Degree Days': PVTt,
            'Stage Sum Degree Days': cumulative_PVTt
        })

        # Check if the current stage threshold is met
        if cumulative_PVTt >= current_stage['Degree_Days'] and current_stage['Degree_Days'] > 0:
            current_stage_index += 1
            if current_stage_index < len(stages):
                current_stage = stages[current_stage_index]
                cumulative_PVTt = 0
                stage_days = 0
            else:
                break

    # Convert results to a DataFrame and save to CSV
    plant_results_df = pd.DataFrame(plant_results)
    output_file = f"True Thermal Time/{plant_id}_True_Thermal_Time.csv"
    plant_results_df.to_csv(output_file, index=False)
    print(f"Saved thermal time data for {plant_id} to {output_file}")
