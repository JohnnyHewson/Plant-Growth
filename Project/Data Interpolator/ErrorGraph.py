import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_data(file_path):
    """Reads data from a CSV file."""
    return pd.read_csv(file_path, parse_dates=['Date'])

def merge_datasets(calculated_data, interpolated_data):
    """Merges the calculated and interpolated datasets on 'Date' and 'Growing Year'."""
    merged_data = pd.merge(
        calculated_data,
        interpolated_data,
        on=['Date', 'Growing Year'],
        suffixes=('_calculated', '_interpolated')
    )
    
    # Calculate Absolute Percentage Error (APE)
    # APE = |(Calculated - Interpolated) / Interpolated| * 100
    # Handle division by zero by setting APE to NaN where interpolated is zero
    merged_data['APE'] = np.where(
        merged_data['Cumulative Thermal Temperature_interpolated'] != 0,
        np.abs((merged_data['Cumulative Thermal Temperature_calculated'] - 
                merged_data['Cumulative Thermal Temperature_interpolated']) / 
               merged_data['Cumulative Thermal Temperature_interpolated']) * 100,
        np.nan
    )
    
    return merged_data

def assign_day_of_growing_season(df):
    """Assigns each date a day number within the growing season (September 1 to August 31)."""
    # Sort the DataFrame by Date
    df = df.sort_values(by='Date')
    
    # Group by Growing Year and assign day numbers
    df['Day_of_Growing_Season'] = df.groupby('Growing Year').cumcount() + 1
    
    return df

def plot_ape(merged_data):
    """Plots the Absolute Percentage Error (APE) for each growing year on the same graph."""
    plt.figure(figsize=(14, 8))
    
    # Get unique growing years
    growing_years = merged_data['Growing Year'].unique()
    
    # Plot APE for each growing year
    for year in growing_years:
        yearly_data = merged_data[merged_data['Growing Year'] == year]
        plt.plot(
            yearly_data['Day_of_Growing_Season'],
            yearly_data['APE'],
            label=year
        )
    
    # Customize the x-axis to show months from September to August
    # Define month start days within the growing season
    # Approximate the day number for each month's start
    month_starts = {
        'Sep': 1,
        'Oct': 31,    # Sep has 30 days
        'Nov': 61,    # Oct has 31 days
        'Dec': 92,    # Nov has 30 days
        'Jan': 122,   # Dec has 31 days
        'Feb': 153,   # Jan has 31 days
        'Mar': 183,   # Feb has 28 or 29 days
        'Apr': 214,   # Mar has 31 days
        'May': 244,   # Apr has 30 days
        'Jun': 275,   # May has 31 days
        'Jul': 305,   # Jun has 30 days
        'Aug': 336    # Jul has 31 days
    }
    
    # Create list of tick positions and labels
    tick_positions = list(month_starts.values())
    tick_labels = list(month_starts.keys())
    
    # Adjust for leap years if necessary (optional)
    # For simplicity, we'll assume non-leap years
    
    plt.xticks(tick_positions, tick_labels)
    
    plt.xlabel('Month')
    plt.ylabel('Absolute Percentage Error (%)')
    plt.title('Absolute Percentage Error in Cumulative Thermal Temperature by Growing Year')
    plt.legend(title='Growing Year', loc='upper left', bbox_to_anchor=(1.05, 1), fontsize='small')
    plt.grid(True)
    plt.tight_layout()
    
    # Show the plot
    plt.show()

def main(calculated_file, interpolated_file):
    # Read the datasets
    calculated_data = read_data(calculated_file)
    interpolated_data = read_data(interpolated_file)
    
    # Merge and calculate APE
    merged_data = merge_datasets(calculated_data, interpolated_data)
    
    # Assign day of growing season
    merged_data = assign_day_of_growing_season(merged_data)
    
    # Plot the Absolute Percentage Error
    plot_ape(merged_data)

# Example usage
if __name__ == "__main__":
    global variable_name
    variable_name = 'ThermalTime'
    year_start = 78  # Starting year in yy format
    year_end = 81    # Ending year in yy format
    calculated_file = f'ThermalTemperature_{year_start}-{year_end}_Calculated.csv'         # Replace with your calculated data file path
    interpolated_file = f"{variable_name}_{year_start}-{year_end}_CumulativeInterpolated.csv"  # Replace with your interpolated data file path
    
    main(calculated_file, interpolated_file)
