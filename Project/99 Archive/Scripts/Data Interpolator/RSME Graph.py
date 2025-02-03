import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def calculate_rmse(series1, series2):
    """Calculate RMSE for two series."""
    return np.sqrt(np.mean((series1 - series2) ** 2))

def main(csv_file1, csv_file2):
    # Load only the required columns and specify the date format
    df1 = pd.read_csv(csv_file1, usecols=['Date', 'Cumulative Thermal Time', 'Growing Year'], parse_dates=['Date'], dayfirst=True)
    df2 = pd.read_csv(csv_file2, usecols=['Date', 'Cumulative Thermal Time', 'Growing Year'], parse_dates=['Date'], dayfirst=True)
    
    # Clean column names (strip spaces, if any)
    df1.columns = df1.columns.str.strip()
    df2.columns = df2.columns.str.strip()
    
    # Merge datasets on 'Date' to ensure alignment
    merged = pd.merge(df1, df2, on='Date', suffixes=('_1', '_2'))
    
    # Clean the column names in the merged dataframe
    merged.columns = merged.columns.str.strip()
    
    # Convert 'Date' column to datetime format (with the correct format for day first)
    merged['Date'] = pd.to_datetime(merged['Date'], format='%d/%m/%Y')  # Specify the correct format
    
    # Extract month and year from the 'Date' column
    merged['Month'] = merged['Date'].dt.month
    merged['Year'] = merged['Date'].dt.year
    
    # Use 'Growing Year_1' to group and calculate RMSE
    rmse_by_month = []
    
    for (year, month), group in merged.groupby(['Growing Year_1', 'Month']):  # Use 'Growing Year_1'
        rmse = calculate_rmse(group['Cumulative Thermal Time_1'], group['Cumulative Thermal Time_2'])
        rmse_by_month.append([year, month, rmse])
    
    # Create a DataFrame with RMSE values per month for each growing year
    rmse_df = pd.DataFrame(rmse_by_month, columns=['Growing Year', 'Month', 'RMSE'])
    
    # Set up month names for plotting
    month_labels = ['S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A']
    
    # Plot RMSE for each growing year
    plt.figure(figsize=(12, 8))
    
    for year, group in rmse_df.groupby('Growing Year'):
        plt.plot(group['Month'], group['RMSE'], label=f'RMSE - {year}')
    
    # Customize the x-axis with month labels
    plt.xticks(range(1, 13), month_labels)  # Use months 1-12 (Sept-Aug)
    plt.xlabel('Month')
    plt.ylabel('RMSE')
    plt.title('RMSE per Month by Growing Year')
    plt.legend(title='Growing Year', loc='upper left', bbox_to_anchor=(1.05, 1), fontsize='small')
    plt.grid(True)
    plt.tight_layout()

    # Show the plot
    plt.show()

# Replace 'file1.csv' and 'file2.csv' with your actual file paths
csv_file1 = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Data Interpolator\\ThermalTime_78-81_CumulativeInterpolated.csv'))
csv_file2 = 'ThermalTime_78-81_Calculated.csv'

main(csv_file1, csv_file2)

