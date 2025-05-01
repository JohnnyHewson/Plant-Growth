import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime

# Load the Interpolated Thermal Time data from the paper
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), r'../..'))
paper_thermal_data_path = os.path.join(project_path, 'Data', 'Interpolated', 'ThermalTime_78-81_CumulativeInterpolated.csv')
paper_thermal_data = pd.read_csv(paper_thermal_data_path, encoding='utf-8')

# # Load additional thermal time data from another source
# additional_thermal_data_path = os.path.join(project_path, 'Data', 'Interpolated', 'ThermalTemperature_78-81_Calculated.csv')
# additional_thermal_data = pd.read_csv(additional_thermal_data_path, encoding='utf-8')

Calculated_Tt_path = os.path.join(project_path,'Data','Processed','Thermal Time')
unaffected_thermal_data = {}
affected_thermal_data = {}
for file in os.listdir(Calculated_Tt_path):
    calc_Tt = pd.read_csv(os.path.join(Calculated_Tt_path,file), encoding='utf-8')
    calc_Tt = calc_Tt[['Date', 'Total Degree Days', 'Sum Unaffected Daily Thermal Time','Stage']]
    calc_Tt['Date'] = pd.to_datetime(calc_Tt['Date'])
    calc_Tt['Growing Year'] = f"{calc_Tt.loc[0,'Date'].year}-{calc_Tt.loc[0,'Date'].year+1}"
    calc_Tt['Day_of_Growing_Season'] = calc_Tt.groupby('Growing Year').cumcount() + (calc_Tt.loc[0,'Date'].timetuple().tm_yday - datetime.datetime(calc_Tt.loc[0,'Date'].year,9,1).timetuple().tm_yday)
    unaffected_thermal_data.update({file[:file.rfind('.')]:calc_Tt[['Date','Growing Year','Day_of_Growing_Season','Sum Unaffected Daily Thermal Time']]})
    affected_thermal_data.update({file[:file.rfind('.')]:calc_Tt[['Date','Growing Year','Day_of_Growing_Season','Total Degree Days','Stage']]})

# Convert Date column to datetime format
paper_thermal_data['Date'] = pd.to_datetime(paper_thermal_data['Date'], format='%d/%m/%Y', dayfirst=True, errors='coerce')
#additional_thermal_data['Date'] = pd.to_datetime(additional_thermal_data['Date'], format='%d/%m/%Y', dayfirst=True, errors='coerce')

# Sort datasets by Date
paper_thermal_data = paper_thermal_data.sort_values(by='Date')
#additional_thermal_data = additional_thermal_data.sort_values(by='Date')

# Calculate day of growing season
paper_thermal_data['Day_of_Growing_Season'] = paper_thermal_data.groupby('Growing Year').cumcount() + 1
#additional_thermal_data['Day_of_Growing_Season'] = additional_thermal_data.groupby('Growing Year').cumcount() + 1

# Set up month markers
month_days = [0, 30, 61, 92, 122, 153, 183, 214, 244, 275, 305, 336]
month_labels = ['S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A']

# Create plot
colors = ['red', 'green', 'blue']  # List of colors
year_colors = {year: colors[i % len(colors)] for i,year in enumerate([1979,1980,1981])}

# Stage abbreviation mapping
stage_labels = {
    'Seeding': 'S',
    'Emergence': 'E',
    'Double Ridge': 'DR',
    'Anthesis': 'A',
    'Maturity': 'M'
}

plt.figure(figsize=(12, 6))
for year, group in paper_thermal_data.groupby('Growing Year'):
    plt.plot(group['Day_of_Growing_Season'], group['Cumulative Thermal Time'], label=f'Paper {year} (Unaffected)', color=year_colors[int(year.split('-')[1])])
# for year, group in additional_thermal_data.groupby('Growing Year'):
#     plt.plot(group['Day_of_Growing_Season'], group['Cumulative Thermal Time'], label=f'Calculated {year}')
for plant_ID in unaffected_thermal_data:
    for year, group in unaffected_thermal_data[plant_ID].groupby('Growing Year'):
        plt.plot(group['Day_of_Growing_Season'], group['Sum Unaffected Daily Thermal Time'], linestyle=':', label=f'Unaffected {plant_ID}', color=year_colors[int(plant_ID.removesuffix('Early' if 'Early' in plant_ID else 'Late'))])
# for plant_ID in affected_thermal_data:
#     for year, group in affected_thermal_data[plant_ID].groupby('Growing Year'):
#         group = group.sort_values(by='Day_of_Growing_Season')
#         plt.plot(group['Day_of_Growing_Season'], group['Total Degree Days'], linestyle='--', label=f'Affected {plant_ID}', color=year_colors[plant_ID])

#         # Find stage changes
#         stage_changes = group[group['Stage'].ne(group['Stage'].shift())]
#         for _, row in stage_changes.iterrows():
#             stage_abbr = stage_labels.get(row['Stage'], row['Stage'])  # fallback to original if not mapped
#             plt.scatter(row['Day_of_Growing_Season'], row['Total Degree Days'], color='black', zorder=5, marker='|')
#             plt.text(row['Day_of_Growing_Season'], row['Total Degree Days'] + 10, stage_abbr,
#                      fontsize=8, ha='center', va='bottom')

plt.xticks(month_days, month_labels)
plt.xlabel('Month')
plt.ylabel('Cumulative Thermal Time')
plt.title('Cumulative Thermal Time by Growing Year')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(project_path,'Data/Graphs',f'ThermalTime.png'))
plt.show()