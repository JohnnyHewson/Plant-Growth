import numpy as np
import pandas as pd
import datetime
import os
from math import pi, cos, radians, sin, sqrt, trunc

configdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config.txt'))
with open(configdir, 'r') as config_file:
    stages = []
    for line in config_file:
        if line.__contains__(','):
            if line[:line.index(',')] == 'shoot_production_rate':
                TPr = float(line.strip().split(',')[1])
            if line[:line.index(',')] == 'Latitude':
                Lat = float(line.strip().split(',')[1])

temp_data = pd.read_csv('Temp 1978-1979.csv')
temp_data['Date'] = pd.to_datetime(temp_data['Date'], dayfirst=True)
path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'Plant Data'))

def calc_RoCoDLatE(latitude, date):
    julian_day = date.timetuple().tm_yday
    tilt_of_earth = 23.44
    radians_per_day = 2 * pi / 365

    latitude_rad = radians(latitude)
    declination = -radians(tilt_of_earth) * cos(radians_per_day * (julian_day + 10))
    declination_rate = (radians(tilt_of_earth) * radians_per_day *
                        sin(radians_per_day * (julian_day + 10)))

    numerator = 24 / pi * cos(declination) * sin(latitude_rad)
    denominator = sqrt(1 - (sin(latitude_rad) * sin(declination))**2)
    
    if denominator == 0:
        return 0.0
    else:
        return numerator / denominator * declination_rate

for file in os.listdir(path):
    plant_data = pd.read_csv(f'Plant Data\\{file}')
    plant_data['Date'] = pd.to_datetime(plant_data['Date'], dayfirst=True)
    plant_data = plant_data.merge(temp_data, on='Date', how='outer')
    plant_data = plant_data.dropna().reset_index().drop(columns='index')
    number_of_growths = pd.DataFrame({'#Tillers':[0],'#Leaves':[0]})

    ### Tiller and Leaf Growth Submodel ###
    rate_of_change_of_daylength_at_emergence = 0
    new_tillers = 0
    new_leaves = 0
    successive_leaf_thermal_time = 0
    offset_total_thermal_time = 0
    max_tillers = 0
    #Constants for Tiller Death Proportion
    A=825
    alpha = 1.46
    beta = 2.24

    for index,row in plant_data.iterrows():
        if row['Stage'] == 'Seeding':
            number_of_growths.loc[index] = [0, 0]
        elif row['Stage'] in ['Emergence','Double Ridge']:
            row['Total Degree Days'] -= offset_total_thermal_time
            if rate_of_change_of_daylength_at_emergence == 0:
                rate_of_change_of_daylength_at_emergence = calc_RoCoDLatE(Lat, row['Date'])
                rate_of_leaf_appearance_per_degree_day = 0.025 * rate_of_change_of_daylength_at_emergence + 0.0104
                phylochron_interval = 1/rate_of_leaf_appearance_per_degree_day
            if number_of_growths['#Leaves'].values[-1] < 3:
                if number_of_growths['#Leaves'].values[-1] == 0:
                    successive_leaf_thermal_time = row['Total Degree Days']
                    new_leaves = rate_of_leaf_appearance_per_degree_day * row['Total Degree Days']
                    number_of_growths.loc[index] = [0, new_leaves]
                else:
                    if (row['Total Degree Days'] - successive_leaf_thermal_time) >= phylochron_interval:
                        successive_leaf_thermal_time = row['Total Degree Days']
                        new_leaves += 1
                    number_of_growths.loc[index] = [0, new_leaves]
            else:
                if row['Stage'] != 'Double Ridge':
                    new_tillers += TPr * sum(plant_data['Daily Degree Days'][index-6:index])
                    if (row['Total Degree Days'] - successive_leaf_thermal_time) >= phylochron_interval:
                        successive_leaf_thermal_time = row['Total Degree Days']
                        new_leaves += 1
                    number_of_growths.loc[index] = [new_tillers, new_leaves]
                else:
                    if (row['Total Degree Days'] - successive_leaf_thermal_time) >= phylochron_interval:
                        successive_leaf_thermal_time = row['Total Degree Days']
                        new_leaves += 1
                    if max_tillers == 0:
                        max_tillers = trunc(number_of_growths['#Tillers'].values[-1])+1
                        N_n = pd.DataFrame({'N_n':[i for i in range(1,max_tillers+1)],
                                            'Survival Chance':[1*(number_of_growths['#Tillers'].values[-1] - trunc(number_of_growths['#Tillers'].values[-1])) if i == trunc(number_of_growths['#Tillers'].values[-1])+1 else 1 for i in range(1,max_tillers+1)]})
                        chance = list(int(i) for i in list('1'*max_tillers))
                    for index2,row2 in N_n.iterrows():
                        chance[index2] *= (1 / (1 + (((min(row['Total Degree Days'],600)/600)/((A/row2['N_n'])**alpha)))**beta))
                    number_of_growths.loc[index] = [new_tillers, new_leaves]             

    N_n = pd.concat([N_n, pd.DataFrame({'Multiplier':chance})], axis=1)
    plant_data = pd.concat([plant_data, number_of_growths], axis=1)

    # Leaf data from the paper (Table 1)
    leaf_data = [
        {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "leaf_area": 653},  # Leaf 1
        {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "leaf_area": 653},  # Leaf 2
        {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "leaf_area": 653},  # Leaf 3
        {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "leaf_area": 653},  # Leaf 4
        {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "leaf_area": 653},  # Leaf 5
        {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "leaf_area": 979},  # Leaf 6
        {"lamina_length": 245, "lamina_width": 11, "sheath_length": 120, "leaf_area": 1797},  # Leaf 7
        {"lamina_length": 275, "lamina_width": 13, "sheath_length": 125, "leaf_area": 2322},  # Leaf 8
        {"lamina_length": 310, "lamina_width": 15, "sheath_length": 125, "leaf_area": 3164},  # Leaf 9
        {"lamina_length": 350, "lamina_width": 15, "sheath_length": 135, "leaf_area": 3481},  # Leaf 10
        {"lamina_length": 395, "lamina_width": 15, "sheath_length": 145, "leaf_area": 3988},  # Leaf 11
        {"lamina_length": 445, "lamina_width": 15, "sheath_length": 155, "leaf_area": 4560},  # Leaf 12
    ]
    projected_area_factor = 1e-6  # Convert from mm^2 to m^2
    
    LAI_z = pd.DataFrame({'Level':[],'LAI':[]})
    LAI = 0
    number_of_levels = 0
    if number_of_growths['#Leaves'].values[-1] >= 12:
            number_of_levels = 5
    else:
        for i in range(trunc(number_of_growths['#Leaves'].values[-1])+1):
            if leaf_data[i]['sheath_length'] > leaf_data[i-1]['sheath_length']:
                number_of_levels += 1
    level = 0
    i=0
    while number_of_levels-level>=0:
        LAI += leaf_data[i]['leaf_area']
        if trunc(number_of_growths['#Leaves'].values[-1] - i) == 0:
            print(number_of_growths['#Leaves'].values[-1], i)
        if i > 0:
            if i < 11:
                LAI += leaf_data[i]['leaf_area']
                if leaf_data[i]['sheath_length'] > leaf_data[i-1]['sheath_length']:
                    level += 1
            else:
                LAI += leaf_data[11]['leaf_area']
            LAI_z.loc[number_of_levels-level] = [f'Level {number_of_levels-level}',LAI]
        i+=1

    # for value in LAI_z:
    #     LAI_z[f'{value}'] *= projected_area_factor
    
    print(LAI_z)
    ### Root Growth Submodel ###

