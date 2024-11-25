import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
from datetime import datetime

# Sample data (replace this with your actual data)

df1 = pd.read_csv("")

data = [
    {"start": "2023-09-01", "end": "2023-10-15", "stage": "S", "dataset": "Dataset 1"},
    {"start": "2023-10-16", "end": "2024-01-10", "stage": "E", "dataset": "Dataset 1"},
    {"start": "2024-01-11", "end": "2024-03-15", "stage": "DR", "dataset": "Dataset 1"},
    {"start": "2024-03-16", "end": "2024-06-01", "stage": "A", "dataset": "Dataset 1"},
    {"start": "2024-06-02", "end": "2024-08-31", "stage": "M", "dataset": "Dataset 1"},
    {"start": "2023-09-05", "end": "2023-10-10", "stage": "S", "dataset": "Dataset 2"},
    {"start": "2023-10-11", "end": "2024-01-05", "stage": "E", "dataset": "Dataset 2"},
    {"start": "2024-01-06", "end": "2024-03-10", "stage": "DR", "dataset": "Dataset 2"},
    {"start": "2024-03-11", "end": "2024-06-05", "stage": "A", "dataset": "Dataset 2"},
    {"start": "2024-06-06", "end": "2024-08-31", "stage": "M", "dataset": "Dataset 2"},
]

# Convert data to DataFrame for easier manipulation
df = pd.DataFrame(data)
df['start'] = pd.to_datetime(df['start'])
df['end'] = pd.to_datetime(df['end'])
df['duration'] = (df['end'] - df['start']).dt.days

# Map stages to colors
stage_colors = {"S": "yellow", "E": "green", "DR": "blue", "A": "orange", "M": "red"}

# Create a figure and axis
fig, ax = plt.subplots(figsize=(12, 6))

# Plot bars for each dataset
datasets = df['dataset'].unique()
y_positions = {dataset: i for i, dataset in enumerate(datasets)}

for _, row in df.iterrows():
    start_date = row['start']
    end_date = row['end']
    y_pos = y_positions[row['dataset']]
    color = stage_colors[row['stage']]
    ax.barh(y_pos, (end_date - start_date).days, left=start_date, color=color, edgecolor='black', label=row['stage'])

# Customize the timeline
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax.set_yticks(range(len(datasets)))
ax.set_yticklabels(datasets)
ax.set_title("Growing Stages Timeline")
ax.set_xlabel("Time")
ax.set_ylabel("Datasets")

# Create a legend
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), title="Stages")

# Show the plot
plt.tight_layout()
plt.show()

#test

#test 2

#test 3