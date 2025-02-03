import pandas as pd
from math import cos, pi, sin, asin, acos, tan
import datetime
import subprocess

# Constants
T_base = 0
T_opt = 26
TD_max = 37

# Load temperature data
year_start = 78
year_end = 81
T_min_df = pd.read_csv(f'Rothamsted Daily Min Temps 01.09.19{year_start}-31.08.19{year_end}.csv')  # CSV with 'Date' and 'T_min' columns
T_max_df = pd.read_csv(f'Rothamsted Daily Max Temps 01.09.19{year_start}-31.08.19{year_end}.csv')  # CSV with 'Date' and 'T_max' columns

# Rename columns to match expected names
T_min_df.rename(columns={'date': 'Date', 'Rothamsted': 'T_min'}, inplace=True)
T_max_df.rename(columns={'date': 'Date', 'Rothamsted': 'T_max'}, inplace=True)

# Convert 'Date' column to datetime format with dayfirst=True
T_min_df['Date'] = pd.to_datetime(T_min_df['Date'], dayfirst=True)
T_max_df['Date'] = pd.to_datetime(T_max_df['Date'], dayfirst=True)

# Merge the dataframes on Date
temperature_df = pd.merge(T_min_df, T_max_df, on="Date")

# Determine the growing year for each date (September 1 to August 31)
temperature_df['Growing_Year'] = temperature_df['Date'].apply(lambda d: d.year if d.month >= 9 else d.year - 1)

# Initialize cumulative Thermal time
cumulative_thermal_time = 0
previous_year = None

# Function to calculate daily Thermal time
def calculate_thermal_time(T_min, T_max):
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

# List to store results
results = []

# Calculate Thermal time for each day
for _, row in temperature_df.iterrows():
    date = row['Date']
    growing_year = row['Growing_Year']
    T_min = max(row['T_min'], 0)
    T_max = max(row['T_max'], 0)

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
        "Growing Year": f"{growing_year}-{growing_year + 1}",
        "Thermal Time": daily_thermal_time,
        "Cumulative Thermal Time": cumulative_thermal_time
    })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Save results to a new CSV file
results_df.to_csv(f'ThermalTime_{year_start}-{year_end}_Calculated.csv', index=False)

subprocess.run(["python", "TtGraphing.py"])