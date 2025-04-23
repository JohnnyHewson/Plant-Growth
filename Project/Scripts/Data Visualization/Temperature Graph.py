import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime

# Load data
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\\..'))
thermal_data_path = os.path.join(project_path, 'Data', 'Raw', 'Temperature 1978-1981.csv')
thermal_data = pd.read_csv(thermal_data_path, encoding='utf-8', header=0)

# Convert Date to datetime
thermal_data['Date'] = pd.to_datetime(thermal_data['Date'])

# Assign Growing Year
def assign_growing_year(date):
    return f"{date.year}-{date.year+1}" if date.month >= 9 else f"{date.year-1}-{date.year}"

thermal_data['Growing Year'] = thermal_data['Date'].apply(assign_growing_year)

# Calculate Day of Growing Season
def day_of_growing_season(row):
    start = datetime.datetime(int(row['Growing Year'].split('-')[0]), 9, 1)
    return (row['Date'] - start).days + 1

thermal_data['Day_of_Growing_Season'] = thermal_data.apply(day_of_growing_season, axis=1)

# Set up subplots for each growing year
growing_years = sorted(thermal_data['Growing Year'].unique())
n_years = len(growing_years)
fig, axes = plt.subplots(n_years, 1, figsize=(12, 3.5 * n_years), sharex=True)

# Month ticks
month_days = [0, 30, 61, 92, 122, 153, 183, 214, 244, 275, 305, 336]
month_labels = ['S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A']

# Plot for each year
for i, year in enumerate(growing_years):
    ax = axes[i] if n_years > 1 else axes
    group = thermal_data[thermal_data['Growing Year'] == year]
    
    ax.fill_between(group['Day_of_Growing_Season'], group['Min Temp'], group['Max Temp'], color='lightpink', label='Temp Range')
    ax.plot(group['Day_of_Growing_Season'], group['Mean Temp'], color='blue', label='Mean Temp')
    
    ax.set_ylabel('Temp (Celcius)')
    ax.set_title(f'Temperature for Growing Year {year}')
    ax.grid(True)
    if i == n_years - 1:
        ax.set_xticks(month_days)
        ax.set_xticklabels(month_labels)
        ax.set_xlabel('Month (Starting in September)')
    else:
        ax.set_xticks([])

    ax.legend()

plt.tight_layout()
plt.show()
