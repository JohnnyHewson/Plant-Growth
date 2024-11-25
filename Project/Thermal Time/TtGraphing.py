import pandas as pd
import matplotlib.pyplot as plt

# Load the Thermal Time and Radiation data
thermal_data = pd.read_csv('ThermalTime_78-81_Calculated.csv', encoding='utf-8')
radiation_data = pd.read_csv('Radiation_78-81_CumulativeInterpolated.csv', encoding='utf-8')

# For thermal data
thermal_data['Date'] = pd.to_datetime(thermal_data['Date'], format='%d/%m/%Y', dayfirst=True, errors='coerce')

# For radiation data
radiation_data['Date'] = pd.to_datetime(radiation_data['Date'], format='%d/%m/%Y', dayfirst=True, errors='coerce')

# Sort both datasets by Date just in case they're not in order
thermal_data = thermal_data.sort_values(by='Date')
radiation_data = radiation_data.sort_values(by='Date')

# Calculate day of growing season for each entry in both datasets
thermal_data['Day_of_Growing_Season'] = thermal_data.groupby('Growing Year').cumcount() + 1
radiation_data['Day_of_Growing_Season'] = radiation_data.groupby('Growing Year').cumcount() + 1

# Set up month markers corresponding to September to August
month_days = [0, 30, 61, 92, 122, 153, 183, 214, 244, 275, 305, 336]
month_labels = ['S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A']

# Create the plot with a second y-axis for radiation data
fig, ax1 = plt.subplots(figsize=(12, 8))

# Plot cumulative thermal time for each growing year
for year, group in thermal_data.groupby('Growing Year'):
    ax1.plot(group['Day_of_Growing_Season'], group['Cumulative Thermal Time'], label=f'Thermal Time {year}')

# Customize the x-axis with month labels
ax1.set_xticks(month_days)
ax1.set_xticklabels(month_labels)
ax1.set_xlabel('Month')
ax1.set_ylabel('Cumulative Thermal Time')
ax1.set_title('Cumulative Thermal Time and Radiation by Growing Year')
ax1.grid(True)

# Create a second y-axis for radiation data
ax2 = ax1.twinx()

# Plot scaled radiation for each growing year on the second axis
for year, group in radiation_data.groupby('Growing Year'):
    ax2.plot(group['Day_of_Growing_Season'], group['Cumulative Radiation'], label=f'Radiation {year}', linestyle='--')

# Set the y-axis label for radiation
ax2.set_ylabel('Scaled Cumulative Radiation')

# Adjust the y-axis limits and ticks for radiation to avoid overlap
ax1.set_ylim(3000)
ax1.set_yticks(range(0, 4000, 1000))
ax2.set_ylim(0, 4000)  # Scale radiation axis to avoid overlap
ax2.set_yticks(range(0, 6001, 2000))  # Adjust number of ticks as needed

# Add legends for both plots
ax1.legend(title='Growing Year (Thermal Time)', loc='upper left', bbox_to_anchor=(1.05, 1), fontsize='small')
ax2.legend(title='Growing Year (Radiation)', loc='upper left', bbox_to_anchor=(1.05, 0.87), fontsize='small')

# Add annotation for the max value of the radiation on the y-axis
radiation_max = radiation_data['Cumulative Radiation'].max()
ax2.annotate(f'Max: {radiation_max:.2f}', xy=(0.95, 0.6), xycoords='axes fraction', 
             ha='center', va='center', fontsize=10, color='black', 
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

# Add annotation for the max value of the thermal time on the y-axis
thermal_max = thermal_data['Cumulative Thermal Time'].max()
ax2.annotate(f'Max: {thermal_max:.2f}', xy=(0.85, 0.95), xycoords='axes fraction', 
             ha='center', va='center', fontsize=10, color='black', 
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

# Tight layout to ensure everything fits
plt.tight_layout()

# Show the plot
plt.show()