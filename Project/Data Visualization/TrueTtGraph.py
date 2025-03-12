import pandas as pd
import matplotlib.pyplot as plt
import os

# Load the Thermal Time and Radiation data
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\\..'))
thermal_data_path = os.path.join(project_path,'Data','Processed','Thermal Time','1979Early Thermal Time.csv') 
thermal_data = pd.read_csv(thermal_data_path, encoding='utf-8')

# For thermal data
thermal_data['Date'] = pd.to_datetime(thermal_data['Date'])

# # Calculate day of growing season for each entry in both datasets
# thermal_data['Day_of_Growing_Season'] = thermal_data.groupby('Growing Year').cumcount() + 1

# Set up month markers corresponding to September to August
month_days = [0, 30, 61, 92, 122, 153, 183, 214, 244, 275, 305, 336]
month_labels = ['S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A']

# Create the plot with a second y-axis for radiation data
fig, ax1 = plt.subplots(figsize=(12, 8))

# # Plot cumulative thermal time for each growing year
# for year, group in thermal_data.groupby('Growing Year'):
#     ax1.plot(group['Day_of_Growing_Season'], group['Cumulative Thermal Time'], label=f'Thermal Time {year}')
ax1.plot(thermal_data['Total Degree Days'], label=f'Thermal Time')

# Customize the x-axis with month labels
ax1.set_xticks(month_days)
ax1.set_xticklabels(month_labels)
ax1.set_xlabel('Month')
ax1.set_ylabel('Cumulative Thermal Time')
ax1.set_title('Cumulative Thermal Time')
ax1.grid(True)

# # Adjust the y-axis limits and ticks for radiation to avoid overlap
# ax1.set_ylim(3000)
# ax1.set_yticks(range(0, 4000, 1000))

# Add legends for both plots
ax1.legend(title='Growing Year (Thermal Time)', loc='upper left', bbox_to_anchor=(1.05, 1), fontsize='small')


# Add annotation for the max value of the thermal time on the y-axis
thermal_max = thermal_data['Total Degree Days'].max()
ax1.annotate(f'Max: {thermal_max:.2f}', xy=(0.85, 0.95), xycoords='axes fraction', 
             ha='center', va='center', fontsize=10, color='black', 
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

# Tight layout to ensure everything fits
plt.tight_layout()

# Show the plot
plt.show()

