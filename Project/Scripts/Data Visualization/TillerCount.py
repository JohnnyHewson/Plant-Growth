import os
import json
import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np

# Get project path two levels up from this file
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), r'../..'))
paper_data_path = os.path.join(project_path,'Data/Interpolated','TillerNumber.csv')

paper_data_messy = pd.read_csv(paper_data_path, header=0)
paper_data = pd.DataFrame()
pairs = [(paper_data_messy.columns[i],paper_data_messy.columns[i+1]) for i in range(0,len(paper_data_messy.columns),2)]
for col1, col2 in pairs:
    paper_data[f'{col1}'] = list(zip(paper_data_messy[col1], paper_data_messy[col2]))


paper_data = paper_data.drop(0,axis=0)

# Load JSON data
with open(os.path.join(project_path, 'Data/Processed', 'graph_data.json')) as file:
    graph_data = json.load(file)

paperDR = {
    '1979Early':200,
    '1980Early':181,
    '1981Early':185,
    '1979Late':225,
    '1980Late':220,
    '1981Late':230
}

# Set up subplot grid
num_plants = len(graph_data)
cols = 3
rows = (num_plants + cols - 1) // cols  # Ceiling division
fig, axes = plt.subplots(rows, cols, figsize=(15, 8), sharey=True)
axes = axes.flatten()

# Plot each plant's data
for i, (plant_name, plant_data) in enumerate(graph_data.items()):
    tillerList = []
    double_ridge_index = None  # Track when stage switches to 'Double Ridge'

    for index in plant_data:
        stage = plant_data[index]['Stage']
        tillerList.append((int(index)+1,sum(plant_data[index]["Number of shoots per m^2"])))
        if double_ridge_index is None and stage == 'Double Ridge':
            double_ridge_index = int(index) + 1  # +1 since x-values are 1-based

    x_values = [x for x,_ in tillerList]
    y_values = [y for _,y in tillerList]

    ax = axes[i]
    ax.plot(x_values, y_values, label='Calculated')

    coords = paper_data[plant_name].map(lambda val: pd.NA if isinstance(val, tuple) and np.isnan(val[0]) and np.isnan(val[1]) else val).dropna()
    coords = coords.sort_values(key=lambda s: s.map(lambda x: x[0]))
    x_paper = [x for x, _ in coords]
    y_paper = [y for _, y in coords]
    ax.plot(x_paper, y_paper, label='Paper')

    # Add vertical dotted line at Double Ridge stage
    if double_ridge_index is not None:
        ax.axvline(x=double_ridge_index, color='C0', linestyle='dotted', linewidth=1)
        ax.text(double_ridge_index + 3, 1920, 'Calculated DR', verticalalignment='center', fontsize=12, color='black')

    ax.axvline(x=paperDR[plant_name], color='C1', linestyle='dotted', linewidth=1)
    ax.text(paperDR[plant_name] + 3, 1790, 'Paper DR', verticalalignment='center', fontsize=12, color='black')
    ax.set_title(plant_name)
    ax.set_yticks([0, 500, 1000, 1500, 2000])
    ax.set_xticks(range(0, 450, 90))
    ax.tick_params(axis='both', labelsize=12)
    ax.grid(axis='y')
    ax.legend(loc='upper left', fontsize=12)

# Hide unused subplots if any
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])
fig.suptitle("Tiller Count against Day of Growing Year", fontsize=16)
fig.supxlabel("Day of Growing Year", fontsize=14)
fig.supylabel("Number of shoots per m²", x=0.01, fontsize=14)
plt.tight_layout()  # leave space for suptitle
plt.savefig(os.path.join(project_path, 'Data/Graphs', 'TillerCount.png'))
plt.show()
plt.close()

