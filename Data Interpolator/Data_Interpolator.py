import pandas as pd
from datetime import datetime, timedelta

def read_data(file_path):
    data = pd.read_csv(file_path, header=None, names=['date', 'value'])
    data['date'] = pd.to_datetime(data['date'], format='%Y/%m/%d')
    data = data[data['value'] >= 0].reset_index(drop=True)
    return data

def generate_full_date_range(start_date, end_date):
    return pd.date_range(start=start_date, end=end_date, freq='D')

def interpolate_values(data, full_dates):
    full_data = pd.DataFrame({'date': full_dates})
    
    # Check if the first date in data is after September 1st
    if data['date'].iloc[0] > full_dates[0]:
        prior_day = full_dates[0] - timedelta(days=1)
        data = pd.concat([pd.DataFrame({'date': [prior_day], 'value': [0]}), data], ignore_index=True)
        full_data = pd.concat([pd.DataFrame({'date': [prior_day]}), full_data], ignore_index=True)

    # Merge the full date DataFrame with the known data
    merged_data = pd.merge(full_data, data, on='date', how='left')

    # Interpolate values forward and backward
    merged_data['value'] = merged_data['value'].interpolate(method='linear')

    if merged_data['value'][0] == 0:
        merged_data.drop([0], inplace=True)

    return merged_data

def process_files(input_files):
    cumulative_data = pd.DataFrame()

    for file_path in input_files:
        year_data = read_data(file_path)
        start_year = year_data['date'].dt.year.iloc[0]
        start_date = datetime(start_year, 9, 1)
        end_date = datetime(start_year + 1, 8, 31)
        
        full_dates = generate_full_date_range(start_date, end_date)
        interpolated_data = interpolate_values(year_data, full_dates)
        
        interpolated_data['Growing Year'] = f"{start_year}-{start_year + 1}"
        
        cumulative_data = pd.concat([cumulative_data, interpolated_data], ignore_index=True)

    cumulative_data.rename(columns={'date': 'Date', 'value': f'Cumulative {variable_name}'}, inplace=True)
    return cumulative_data[['Date', 'Growing Year', f'Cumulative {variable_name}']]

def write_data(output_file_path, data):
    data.to_csv(output_file_path, index=False)

def main(input_files, output_file):
    # Process each input file and combine the data
    cumulative_data = process_files(input_files)
    
    # Write the result to the output file
    write_data(output_file, cumulative_data)

# Example usage
global variable_name
variable_name = 'Radiation'
year_start = 78  # Starting year in yy format
year_end = 81    # Ending year in yy format

input_files = [f"{variable_name}\\{variable_name}{year}-" + f"{year + 1}Cumulative.csv" for year in range(year_start, year_end)]
output_file = f"{variable_name}_{year_start}-{year_end}_CumulativeInterpolated.csv"

main(input_files, output_file)
