import pandas as pd
import matplotlib.pyplot as plt
import os

# Load the Interpolated Thermal Time data from the paper
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\\..'))
paper_thermal_data_path = os.path.join(project_path, 'Data', 'Interpolated', 'ThermalTime_78-81_CumulativeInterpolated.csv')
paper_thermal_data = pd.read_csv(paper_thermal_data_path, encoding='utf-8')

# Load additional thermal time data from another source
additional_thermal_data_path = os.path.join(project_path, 'Data', 'Interpolated', 'ThermalTemperature_78-81_Calculated.csv')
additional_thermal_data = pd.read_csv(additional_thermal_data_path, encoding='utf-8')

# Convert Date column to datetime format
paper_thermal_data['Date'] = pd.to_datetime(paper_thermal_data['Date'], format='%d/%m/%Y', dayfirst=True, errors='coerce')
additional_thermal_data['Date'] = pd.to_datetime(additional_thermal_data['Date'], format='%d/%m/%Y', dayfirst=True, errors='coerce')

# Sort datasets by Date
paper_thermal_data = paper_thermal_data.sort_values(by='Date')
additional_thermal_data = additional_thermal_data.sort_values(by='Date')

# Calculate day of growing season
paper_thermal_data['Day_of_Growing_Season'] = paper_thermal_data.groupby('Growing Year').cumcount() + 1
additional_thermal_data['Day_of_Growing_Season'] = additional_thermal_data.groupby('Growing Year').cumcount() + 1

# Set up month markers
month_days = [0, 30, 61, 92, 122, 153, 183, 214, 244, 275, 305, 336]
month_labels = ['S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A']

# Create plot
plt.figure(figsize=(12, 6))
for year, group in paper_thermal_data.groupby('Growing Year'):
    plt.plot(group['Day_of_Growing_Season'], group['Cumulative Thermal Time'], label=f'Paper {year}')
for year, group in additional_thermal_data.groupby('Growing Year'):
    plt.plot(group['Day_of_Growing_Season'], group['Cumulative Thermal Time'], linestyle='-.', label=f'Calculated {year}')

plt.xticks(month_days, month_labels)
plt.xlabel('Month')
plt.ylabel('Cumulative Thermal Time')
plt.title('Cumulative Thermal Time by Growing Year')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()