import os
import matplotlib.pyplot as plt
import json

data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), r'../..','Data'))

with open (os.path.join(data_path,'Processed','graph_data.json')) as graph_data_file:
    graph_data = json.load(graph_data_file)

LAIdata = []

for i,PlantID in enumerate(graph_data):
    LAIdata.append((PlantID,[]))
    for index in graph_data[PlantID]:
        LAIdata[i][1].append(graph_data[PlantID][index]['LAI'])

fig, axs = plt.subplots(2, 3)
axs = axs.flatten()
for n,LAIGraph in enumerate(LAIdata):
    axs[n].plot([x for x,_ in LAIGraph[1]],[y for _,y in LAIGraph[1]])
    axs[n].set_xticks([0,90,180,270,360])
    axs[n].set_yticks([0,2.5,5,7.5,10])
    axs[n].set_title(f'Plant {LAIGraph[0]}')
    axs[n].set_xlabel('Days Since Planting')
    axs[n].set_ylabel('Peak LAI')
    axs[n].grid()

fig.suptitle("Peak LAI vs. Days Since Planting", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(os.path.join(data_path,'Graphs','Peak_LAI.png'))
plt.show()
