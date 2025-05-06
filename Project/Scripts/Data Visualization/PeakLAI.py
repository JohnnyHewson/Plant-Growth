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

fig, axs = plt.subplots(2, 3,sharex=True,sharey=True)
axs = axs.flatten()
for n,(PlantID,LAIGraph) in enumerate(LAIdata):
    axs[n].plot([x for x,_ in LAIGraph],[y for _,y in LAIGraph])
    axs[n].set_xticks([0,90,180,270,360])
    axs[n].set_yticks([0,2.5,5,7.5,10])
    axs[n].tick_params('both',labelsize=12)
    axs[n].set_title(PlantID,fontsize=12)
    axs[n].grid()

fig.suptitle("LAI vs. Days Since Planting", fontsize=16)
fig.supxlabel("Day of Growing Year", fontsize=14)
fig.supylabel("LAI", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(data_path,'Graphs','Peak_LAI.png'))
plt.show()
