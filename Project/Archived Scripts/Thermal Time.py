import pandas as pd
from math import cos, pi, sin, asin, acos, tan
import datetime
import os

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
while project_path[-7:] != "Project":
    project_path = os.path.abspath(os.path.join(project_path, '..'))

# Constants
T_base = 0
T_opt = 26
TD_max = 37

# Load temperature data
year_start = 78
year_end = 81
temperature_df = pd.read_csv(os.path.join(project_path, 'Data', 'Raw','Temperature 1978-1981.csv'))
temperature_df['Date'] = pd.to_datetime(temperature_df['Date'])
for i in range(len(temperature_df)):
    temperature_df.loc[i,'Growing Year'] = temperature_df.loc[i,'Date'].year if temperature_df.loc[i,'Date'].month >= 9 else temperature_df.loc[i,'Date'].year - 1

# Initialize cumulative Thermal time
cumulative_thermal_time = 0
previous_year = None

# Function to calculate daily Thermal time
def calculate_thermal_time(T_min, T_max):
    T_t = 0
    T_min = max(T_min,0)
    T_max = max(T_max,0)
    for r in range(1, 9):
        f_r = (1 / 2) * (1 + cos((90 / 8) * (2 * r - 1) * pi / 180))
        T_H = max(T_min + f_r * (T_max - T_min), 0)  # Degree Celsius
        if T_H < T_opt:
            T_t += T_H - T_base
        elif T_H == T_opt:
            T_t += T_opt - T_base
        elif T_opt < TD_max:
            T_t += (T_opt - T_base) * (TD_max - T_H) / (TD_max - T_opt)
        else:
            T_t += 0
    return max((1 / 8) * T_t, 0)  # Degree Celsius Days

# List to store results
results = []

# Calculate Thermal time for each day
for _, row in temperature_df.iterrows():
    date = row['Date']
    growing_year = row['Growing Year']
    T_min = max(row['Min Temp'], 0)
    T_max = max(row['Max Temp'], 0)

    # Reset cumulative Thermal time at the start of each new growing year
    if previous_year is not None and growing_year != previous_year:
        cumulative_thermal_time = 0
    previous_year = growing_year
    
    # Calculate the day's Thermal time
    daily_thermal_time = calculate_thermal_time(T_min, T_max)
    
    # Update cumulative Thermal time
    cumulative_thermal_time += daily_thermal_time
    
    # Append results
    results.append({
        "Date": date,
        "Growing Year": f'{growing_year}-{growing_year+1}',
        "Thermal Time": daily_thermal_time,
        "Cumulative Thermal Time": cumulative_thermal_time
    })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Save results to a new CSV file
save_directory = os.path.join(project_path,'Archive','Data','Base 1 Thermal Time',f'ThermalTime_{year_start}-{year_end}_Calculated.csv')
results_df.to_csv(save_directory, index=False)
print("Complete")