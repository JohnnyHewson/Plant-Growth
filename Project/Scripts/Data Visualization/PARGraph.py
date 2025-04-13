import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Load PAR data (replace 'hourly_par_data.csv' with your data file)
# Assumes CSV format with 24 columns (0-23 hours) and 365 rows (days)
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\\..'))
par_path = os.path.join(project_path,'Data','Processed','Average Conditions','Average Hourly Par.csv')
par_data = pd.read_csv(par_path, index_col=0)

# Conversion factor: 1 W/m^2 for 1 hour = 3.6 kJ/m^2 = 0.0036 MJ/m^2
W_TO_MJ_CONVERSION = 3600 / 1e6  # (seconds/hour / J/MJ)

# Clean column names (convert to integers)
par_data.columns = par_data.columns.astype(int)
par_data = par_data * W_TO_MJ_CONVERSION

# Create coordinate grids
X = np.arange(25)  # Hours 0-24 (needs +1 for pcolormesh boundaries)
Y = np.arange(par_data.shape[0] + 1)  # Days boundaries

# Create Daily Hourly PAR plot (heatmap)
plt.figure(figsize=(12, 8))
plt.pcolormesh(X, Y, par_data.values,
               shading='flat',
               cmap='viridis')
plt.colorbar(label='PAR (MJ/m^2)')
plt.xlabel('Hour of the Day')
plt.ylabel('Day of the Year')
plt.title('Average Hourly PAR')
plt.xticks(np.arange(0, 25, 3))
plt.yticks(np.arange(0, 367, 50))
plt.gca().invert_yaxis()
plt.tight_layout()

# Create Daily Cumulative PAR plot
daily_cumulative = par_data.sum(axis=1)

plt.figure(figsize=(12, 6))
plt.plot(np.arange(1, len(par_data) + 1), daily_cumulative)
plt.xlabel('Day of the Year')
plt.ylabel('Cumulative PAR (MJ/m^2)')
plt.title('Average Daily PAR')
plt.grid(True)
plt.xlim(1, len(par_data))
plt.xticks(np.arange(0, 367, 50))
plt.tight_layout()

# Calculate and plot running total in MJ/m^2
daily_sums = par_data.sum(axis=1)
running_total = daily_sums.cumsum()

plt.figure(figsize=(12, 6))
plt.plot(np.arange(1, len(par_data) + 1), running_total, color='darkgreen')
plt.xlabel('Day of the Year')
plt.ylabel('Cumulative PAR (MJ/m^2)')
plt.title('Average Annual PAR')
plt.grid(True)
plt.xlim(1, len(par_data))
plt.xticks(np.arange(0, len(par_data)+1, 50))
plt.tight_layout()

plt.show()