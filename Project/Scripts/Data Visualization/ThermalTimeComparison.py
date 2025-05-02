import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime

# Load the Interpolated Thermal Time data from the paper
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), r'../..'))
paper_thermal_data_path = os.path.join(project_path, 'Data', 'Interpolated', 'ThermalTime_78-81_CumulativeInterpolated.csv')
paper_thermal_data = pd.read_csv(paper_thermal_data_path, encoding='utf-8')

# Convert Date column to datetime format
paper_thermal_data['Date'] = pd.to_datetime(paper_thermal_data['Date'], format='%d/%m/%Y', errors='coerce')

# Sort datasets by Date
paper_thermal_data = paper_thermal_data.sort_values(by='Date')


# Calculate day of growing season
paper_thermal_data['Day_of_Growing_Season'] = paper_thermal_data.groupby('Growing Year').cumcount() + 1


Calculated_Tt_path = os.path.join(project_path,'Data','Processed','Thermal Time')
unaffected_thermal_data = {}
affected_thermal_data = {}
for file in os.listdir(Calculated_Tt_path):
    calc_Tt = pd.read_csv(os.path.join(Calculated_Tt_path,file), encoding='utf-8')
    calc_Tt = calc_Tt[['Date', 'Total Degree Days', 'Sum Unaffected Daily Thermal Time','Stage']]
    calc_Tt['Date'] = pd.to_datetime(calc_Tt['Date'])
    calc_Tt['Sum Unaffected Daily Thermal Time'] = calc_Tt['Sum Unaffected Daily Thermal Time'].add(paper_thermal_data.loc[calc_Tt.loc[0,'Date'] == paper_thermal_data['Date'],'Cumulative Thermal Time'].values[0] - calc_Tt.loc[0,'Sum Unaffected Daily Thermal Time'])
    calc_Tt['Growing Year'] = f"{calc_Tt.loc[0,'Date'].year}-{calc_Tt.loc[0,'Date'].year+1}"
    calc_Tt['Day_of_Growing_Season'] = calc_Tt.groupby('Growing Year').cumcount() + (calc_Tt.loc[0,'Date'].timetuple().tm_yday - datetime.datetime(calc_Tt.loc[0,'Date'].year,9,1).timetuple().tm_yday)
    unaffected_thermal_data.update({file[:file.rfind('.')]:calc_Tt[['Date','Growing Year','Day_of_Growing_Season','Sum Unaffected Daily Thermal Time']]})
    affected_thermal_data.update({file[:file.rfind('.')]:calc_Tt[['Date','Growing Year','Day_of_Growing_Season','Total Degree Days','Stage']]})


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

# Create plot per growing year
for year in paper_thermal_data['Growing Year'].unique():
    plt.figure(figsize=(12, 6))
    
    # Plot paper data for this year
    group = paper_thermal_data[paper_thermal_data['Growing Year'] == year]
    year_color = year_colors[int(year.split('-')[1])]
    plt.plot(group['Day_of_Growing_Season'], group['Cumulative Thermal Time'],
             label=f'Paper {year} (Unaffected)', color=year_color)
    
    # Plot all relevant plant data for this year
    for plant_ID in unaffected_thermal_data:
        plant_data = unaffected_thermal_data[plant_ID]
        if year in plant_data['Growing Year'].values:
            group = plant_data[plant_data['Growing Year'] == year]
            plant_year = int(plant_ID.removesuffix('Early' if 'Early' in plant_ID else 'Late'))
            plant_color = year_colors.get(plant_year, 'gray')
            plt.plot(group['Day_of_Growing_Season'], group['Sum Unaffected Daily Thermal Time'],
                     linestyle=':', label=f'Unaffected {plant_ID}', color=plant_color)

    # Finalize the plot
    plt.xticks(month_days, month_labels)
    plt.xlabel('Month')
    plt.ylabel('Cumulative Thermal Time')
    plt.title(f'Cumulative Thermal Time - Growing Year {year}')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(project_path, 'Data/Graphs', f'ThermalTime_{year}.png'))
    plt.show()
    plt.close()

# Set up output directory for difference plots
difference_plot_dir = os.path.join(project_path, 'Data/Graphs/DifferencePlots')
os.makedirs(difference_plot_dir, exist_ok=True)
# Set up output directory for percent error plots
percent_error_plot_dir = os.path.join(project_path, 'Data/Graphs/PercentErrorPlots')
os.makedirs(percent_error_plot_dir, exist_ok=True)
# Output directory
stacked_dir = os.path.join(project_path, 'Data/Graphs', 'StackedPlots')
os.makedirs(stacked_dir, exist_ok=True)

for year in paper_thermal_data['Growing Year'].unique():
    paper_group = paper_thermal_data[paper_thermal_data['Growing Year'] == year][['Date', 'Cumulative Thermal Time']]

    for plant_ID, plant_df in unaffected_thermal_data.items():
        if year in plant_df['Growing Year'].values:
            sim_group = plant_df[plant_df['Growing Year'] == year][['Date', 'Sum Unaffected Daily Thermal Time']]
            merged = pd.merge(paper_group, sim_group, on='Date', how='inner')
            if merged.empty:
                continue

            merged['Difference'] = merged['Sum Unaffected Daily Thermal Time'] - merged['Cumulative Thermal Time']

            # Plot and save
            plt.figure(figsize=(10, 4))
            plt.plot(merged['Date'], merged['Difference'], label=f'Difference {plant_ID}')
            plt.axhline(0, color='gray', linestyle='--')
            plt.title(f'Difference in Thermal Time - {plant_ID} ({year})')
            plt.ylabel('Simulated - Observed')
            plt.xlabel('Date')
            plt.legend()
            plt.grid(True)
            plt.tight_layout()

            # File-safe name
            safe_year = year.replace('/', '-').replace(':', '-')
            filename = f'Diff_{plant_ID}_{safe_year}.png'
            filepath = os.path.join(difference_plot_dir, filename)
            plt.savefig(filepath)
            plt.close()  # Close the plot to free memory

for year in paper_thermal_data['Growing Year'].unique():
    paper_group = paper_thermal_data[paper_thermal_data['Growing Year'] == year][['Date', 'Cumulative Thermal Time', 'Day_of_Growing_Season']]

    plant_ids = [pid for pid, df in unaffected_thermal_data.items() if year in df['Growing Year'].values]
    if not plant_ids:
        continue

    for plant_ID in plant_ids:
        sim_group = unaffected_thermal_data[plant_ID]
        sim_year = sim_group[sim_group['Growing Year'] == year][['Date', 'Sum Unaffected Daily Thermal Time', 'Day_of_Growing_Season']]
        merged = pd.merge(paper_group, sim_year, on='Date', how='inner')

        if merged.empty or (merged['Cumulative Thermal Time'] == 0).all():
            continue

        # Calculate differences and errors
        merged['Difference'] = merged['Sum Unaffected Daily Thermal Time'] - merged['Cumulative Thermal Time']
        merged = merged[merged['Cumulative Thermal Time'] != 0]
        merged['Percent Error'] = (merged['Difference'] / merged['Cumulative Thermal Time']) * 100

        # Set up figure
        fig, axs = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
        fig.suptitle(f'Thermal Time Comparison - {plant_ID} ({year})', fontsize=16)

        # 1. Observed vs Simulated
        # Recalculate Day_of_Growing_Season after merge
        start_date = merged['Date'].min()
        merged['Day_of_Growing_Season'] = (merged['Date'] - start_date).dt.days + 1
        axs[0].plot(merged['Day_of_Growing_Season'], merged['Cumulative Thermal Time'], label='Observed', color='black')
        axs[0].plot(merged['Day_of_Growing_Season'], merged['Sum Unaffected Daily Thermal Time'], label='Simulated', linestyle='--', color='blue')
        axs[0].set_ylabel('Cumulative TT')
        axs[0].legend()
        axs[0].grid(True)
        axs[0].set_title('Observed vs Simulated Thermal Time')

        # 2. Difference plot
        axs[1].plot(merged['Day_of_Growing_Season'], merged['Difference'], color='orange')
        axs[1].axhline(0, color='gray', linestyle='--')
        axs[1].set_ylabel('Difference')
        axs[1].grid(True)
        axs[1].set_title('Difference (Simulated - Observed)')

        # 3. Percent error
        axs[2].plot(merged['Day_of_Growing_Season'], merged['Percent Error'], color='red')
        axs[2].axhline(0, color='gray', linestyle='--')
        axs[2].set_ylabel('% Error')
        axs[2].set_xlabel('Day of Growing Season')
        axs[2].grid(True)
        axs[2].set_title('Percent Error')

        # Save and close
        filename = f'StackedComparison_{plant_ID}_{year}.png'
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(os.path.join(stacked_dir, filename))
        plt.close()


