import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os

data = []
# Sample data (replace this with your actual data)
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), r'../..'))
path = os.path.join(project_path, 'Data', 'Processed', 'Thermal Time')
for file in os.listdir(path):
    plant_ID = file.split(' ')[0]
    plant_data = pd.read_csv(os.path.join(path, file))
    plant_data['Date'] = pd.to_datetime(plant_data['Date'])
    group = plant_data.groupby("Stage", sort=False)["Date"]
    agg = group.aggregate(['first', 'last', 'size'])
    for stage in agg.index:
        data.append({
            "start": agg.loc[stage, "first"],
            "end": agg.loc[stage, "last"],
            "duration": agg.loc[stage, "size"],
            "stage": stage,
            "dataset": plant_ID
        })

# Convert data to DataFrame for easier manipulation
df = pd.DataFrame(data)
df['start'] = pd.to_datetime(df['start'])
df['end'] = pd.to_datetime(df['end'])

# Normalize dates to same synthetic year: Sept 1 -> Aug 31 of next year
def normalize_date(date):
    # Map Sept-Dec to 2023, Jan-Aug to 2024
    if date.month >= 9:
        return date.replace(year=2023)
    else:
        return date.replace(year=2024)

df['norm_start'] = df['start'].apply(normalize_date)
df['norm_end'] = df['end'].apply(normalize_date)

# Map stages to colors
stage_colors = {
    "Seeding": "yellow",
    "Emergence": "green",
    "Double Ridge": "blue",
    "Anthesis": "orange",
    "Maturity": "red"
}

# Create a figure and axis
fig, ax = plt.subplots(figsize=(12, 6))

# Assign Y positions
datasets = df['dataset'].unique()
y_positions = {dataset: i for i, dataset in enumerate(datasets)}

# Plot bars
for _, row in df.iterrows():
    y_pos = y_positions[row['dataset']]
    color = stage_colors.get(row['stage'], 'gray')
    duration = (row['norm_end'] - row['norm_start']).days
    ax.barh(
        y_pos,
        duration,
        left=row['norm_start'],
        color=color,
        edgecolor='black',
        label=row['stage']
    )

# Format x-axis: September to August
ax.set_xlim(datetime(2023, 9, 1), datetime(2024, 8, 31))
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))

# Y-axis setup
ax.set_yticks(range(len(datasets)))
ax.set_yticklabels(str(plant).removesuffix('.csv') for plant in datasets)

# Labels and title
ax.set_title("Growing Stages Timeline")
ax.set_xlabel("Month")
ax.set_ylabel("Plant ID")

# Remove duplicate legend labels
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), title="Stages")


plt.tight_layout()
plt.savefig(os.path.join(project_path,'Data/Graphs',f'Stage_Timings.png'))
plt.show()
