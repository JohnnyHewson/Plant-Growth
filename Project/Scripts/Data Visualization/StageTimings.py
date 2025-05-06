import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import json

# Load data
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), r'../..'))
paper_results_path = os.path.join(project_path, 'Data/Interpolated', 'Observed and Simulated Stage Timings.json')
with open(paper_results_path) as paper:
    StageTimings = json.load(paper)

for plant in StageTimings:
    StageTimings[plant].update({'Seeding Date': datetime.strptime(StageTimings[plant]['Seeding Date'], '%Y/%m/%d')})

calc_St_path = os.path.join(project_path, 'Data/Processed/Thermal Time')
for file in os.listdir(calc_St_path):
    plant_ID = file.removesuffix('.csv')
    plant_data = pd.read_csv(os.path.join(calc_St_path, file))
    plant_data['Date'] = pd.to_datetime(plant_data['Date'])
    group = plant_data.groupby("Stage", sort=False)["Date"]
    agg = group.aggregate(['first', 'last', 'size'])
    StageTimings[plant_ID]['Calculated'] = {
        'Seeding': int(agg.loc['Seeding', 'size']),
        'Emergence': int(agg.loc['Emergence', 'size']),
        'Double Ridge': int(agg.loc['Double Ridge', 'size']),
        'Anthesis': int(agg.loc['Anthesis', 'size'])
    }

# Setup plot
fig, ax = plt.subplots(figsize=(10, 9))  # 50% taller

# Stage markers and colors
stages = ['Seeding', 'Emergence', 'Double Ridge', 'Anthesis']
colors = {'Observed': 'black', 'Paper': 'blue', 'Calculated': 'red'}
linestyle = {'Observed': '-', 'Paper': '--', 'Calculated': '-.'}

# Timeline base
y_offset = (len(StageTimings.items()) - 5) * 1.5  # Spread initial offset
yticks = []
ytick_labels = []

for entry_name, entry in StageTimings.items():
    seed_date = datetime(1 if entry["Seeding Date"].month >= 9 else 2, entry["Seeding Date"].month, entry["Seeding Date"].day)

    for method in ['Observed', 'Paper', 'Calculated']:
        ax.plot(seed_date, y_offset, marker='|', markersize=8, color='black')
        for stage in stages:
            if stage in ['Seeding', 'Emergence', 'Anthesis']:
                ax.annotate(text='S' if stage == 'Seeding'
                            else 'E' if stage == 'Emergence'
                            else 'A',
                            xy=(seed_date if stage == 'Seeding' else date, y_offset),
                            xytext=(seed_date - timedelta(days=2) if stage == 'Seeding' else date - timedelta(days=2), y_offset + 0.35),
                            fontsize=12)
            else:
                ax.annotate(text='DR', xy=(seed_date if stage == 'Seeding' else date, y_offset),
                            xytext=(date - timedelta(days=5), y_offset + 0.35), fontsize=12)
            day_offset = entry[method][stage]
            if stage == 'Seeding':
                date = seed_date + timedelta(days=day_offset)
            else:
                date += timedelta(days=day_offset)
            ax.plot(date, y_offset, marker='|', markersize=8, color='black')

            # Horizontal reference line
            ax.hlines(y_offset, seed_date, date, colors='black', linestyle=linestyle[method])
        ax.annotate(text='M', xy=(date, y_offset), xytext=(date - timedelta(days=3), y_offset + 0.35), fontsize=12)
        yticks.append(y_offset)
        ytick_labels.append(entry_name + ' ' + method)
        if method == 'Calculated' and entry_name != '1981Late':
            ax.hlines(y_offset - 0.45, datetime(1, 9, 1), datetime(2, 9, 1), colors='black')
        y_offset -= 1.5  # More space between rows

# Formatting
ax.set_yticks(yticks)
ax.set_yticklabels([label.split(' ')[-1] for label in ytick_labels], fontsize=14)
ax.tick_params('x', labelsize=14)

# Group labels
for y, label in zip(yticks, ytick_labels):
    if 'Paper' in label:
        ax.text(-0.16, y, label[:5], va='center', ha='right',
                fontsize=14, transform=ax.get_yaxis_transform(), fontweight='bold')

ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.set_xlim(datetime(1, 9, 1), datetime(2, 9, 1))
ax.set_xlabel("Month", fontsize=14)
ax.set_title("Phenological Stage Comparison", fontsize=16)
ax.grid(True, axis='x', linestyle='-')
plt.subplots_adjust(left=0.1)
plt.tight_layout()
plt.savefig(os.path.join(project_path, 'Data/Graphs', 'StageTimings.png'))
plt.show()
