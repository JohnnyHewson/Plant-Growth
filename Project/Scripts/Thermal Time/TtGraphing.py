import pandas as pd
import matplotlib.pyplot as plt
import os

#Load the Interpolated Thermal Time and Radiation data from the paper
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\\..'))
paper_thermal_data_path = os.path.join(project_path,'Data','Interpolated','ThermalTime_78-81_CumulativeInterpolated.csv') 
paper_thermal_data = pd.read_csv(paper_thermal_data_path, encoding='utf-8')
paper_radiation_data_path = os.path.join(project_path,'Data','Interpolated','Radiation_78-81_CumulativeInterpolated.csv')
paper_radiation_data = pd.read_csv(paper_radiation_data_path, encoding='utf-8')

#Load the Calculatted Thermal Time and Radiation data from external databases
thermal_data_path = os.path.join(project_path,'Data','Processed','Thermal Time','1979Early Thermal Time.csv') 
thermal_data = pd.read_csv(thermal_data_path, encoding='utf-8')
paper_radiation_data_path = os.path.join(project_path,'Data','Interpolated','Radiation_78-81_CumulativeInterpolated.csv')
paper_radiation_data = pd.read_csv(paper_radiation_data_path, encoding='utf-8')

#For thermal data
paper_thermal_data['Date'] = pd.to_datetime(paper_thermal_data['Date'], format='%d/%m/%Y', dayfirst=True, errors='coerce')

#For radiation data
paper_radiation_data['Date'] = pd.to_datetime(paper_radiation_data['Date'], format='%d/%m/%Y', dayfirst=True, errors='coerce')

#Sort both datasets by Date just in case they're not in order
paper_thermal_data = paper_thermal_data.sort_values(by='Date')
paper_radiation_data = paper_radiation_data.sort_values(by='Date')

#Calculate day of growing season for each entry in both datasets
paper_thermal_data['Day_of_Growing_Season'] = paper_thermal_data.groupby('Growing Year').cumcount() + 1
paper_radiation_data['Day_of_Growing_Season'] = paper_radiation_data.groupby('Growing Year').cumcount() + 1

#Set up month markers corresponding to September to August
month_days = [0, 30, 61, 92, 122, 153, 183, 214, 244, 275, 305, 336]
month_labels = ['S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A']

#Create the plot with a second y-axis for radiation data
fig, ax1 = plt.subplots(figsize=(12, 8))

#Plot cumulative thermal time for each growing year
for year, group in paper_thermal_data.groupby('Growing Year'):
    ax1.plot(group['Day_of_Growing_Season'], group['Cumulative Thermal Time'], label=f'Thermal Time {year}')

#Customize the x-axis with month labels
ax1.set_xticks(month_days)
ax1.set_xticklabels(month_labels)
ax1.set_xlabel('Month')
ax1.set_ylabel('Cumulative Thermal Time')
ax1.set_title('Cumulative Thermal Time and Radiation by Growing Year')
ax1.grid(True)

#Create a second y-axis for radiation data
ax2 = ax1.twinx()

#Plot scaled radiation for each growing year on the second axis
for year, group in paper_radiation_data.groupby('Growing Year'):
    ax2.plot(group['Day_of_Growing_Season'], group['Cumulative Radiation'], label=f'Radiation {year}', linestyle='--')

#Set the y-axis label for radiation
ax2.set_ylabel('Scaled Cumulative Radiation')

#Adjust the y-axis limits and ticks for radiation to avoid overlap
ax1.set_ylim(3000)
ax1.set_yticks(range(0, 4000, 1000))
ax2.set_ylim(0, 4000)  # Scale radiation axis to avoid overlap
ax2.set_yticks(range(0, 6001, 2000))  # Adjust number of ticks as needed

#Add legends for both plots
ax1.legend(title='Growing Year (Thermal Time)', loc='upper left', bbox_to_anchor=(1.05, 1), fontsize='small')
ax2.legend(title='Growing Year (Radiation)', loc='upper left', bbox_to_anchor=(1.05, 0.87), fontsize='small')

#Add annotation for the max value of the radiation on the y-axis
radiation_max = paper_radiation_data['Cumulative Radiation'].max()
ax2.annotate(f'Max: {radiation_max:.2f}', xy=(0.95, 0.6), xycoords='axes fraction', 
             ha='center', va='center', fontsize=10, color='black', 
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

#Add annotation for the max value of the thermal time on the y-axis
thermal_max = paper_thermal_data['Cumulative Thermal Time'].max()
ax2.annotate(f'Max: {thermal_max:.2f}', xy=(0.85, 0.95), xycoords='axes fraction', 
             ha='center', va='center', fontsize=10, color='black', 
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

#Tight layout to ensure everything fits
plt.tight_layout()

#Show the plot
plt.show()
