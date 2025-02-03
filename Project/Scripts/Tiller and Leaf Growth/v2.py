from logging import root
import numpy as np
import pandas as pd
import datetime
import os
from math import pi, cos, radians, sin, sqrt, trunc, exp

configdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config.txt'))
with open(configdir, 'r') as config_file:
    stages = []
    for line in config_file:
        if line.__contains__(','):
            if line[:line.index(',')] == 'shoot_production_rate':
                TPr = float(line.strip().split(',')[1])
            if line[:line.index(',')] == 'Latitude':
                Lat = float(line.strip().split(',')[1])
            if line[:line.index(',')] == 'extinction_coefficient':
                k = float(line.strip().split(',')[1])
            if line[:line.index(',')] == 'leaf_transmission_coefficient':
                m = float(line.strip().split(',')[1])

# Leaf data from the paper (Table 1)
leaf_data = [
    {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "max_leaf_area": 653},  # Leaf 1
    {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "max_leaf_area": 653},  # Leaf 2
    {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "max_leaf_area": 653},  # Leaf 3
    {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "max_leaf_area": 653},  # Leaf 4
    {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "max_leaf_area": 653},  # Leaf 5
    {"lamina_length": 125, "lamina_width": 6, "sheath_length": 0, "max_leaf_area": 979},  # Leaf 6
    {"lamina_length": 245, "lamina_width": 11, "sheath_length": 120, "max_leaf_area": 1797},  # Leaf 7
    {"lamina_length": 275, "lamina_width": 13, "sheath_length": 125, "max_leaf_area": 2322},  # Leaf 8
    {"lamina_length": 310, "lamina_width": 15, "sheath_length": 125, "max_leaf_area": 3164},  # Leaf 9
    {"lamina_length": 350, "lamina_width": 15, "sheath_length": 135, "max_leaf_area": 3481},  # Leaf 10
    {"lamina_length": 395, "lamina_width": 15, "sheath_length": 145, "max_leaf_area": 3988},  # Leaf 11
    {"lamina_length": 445, "lamina_width": 15, "sheath_length": 155, "max_leaf_area": 4560},  # Leaf 12
]
projected_area_factor = 1e-6  # Convert from mm^2 to m^2

path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'Plant Data'))

temp_data = pd.read_csv('Temp 1978-1979.csv')
temp_data['Date'] = pd.to_datetime(temp_data['Date'], dayfirst=True)
temp_data = temp_data.rename(columns={'Min_Temp':'Min Temp','Max_Temp':'Max Temp','Mean_Temp':'Mean Temp'})
PAR_data = pd.read_csv('Average Hourly PAR.csv')

# Function to calculate daily thermal time (degree days)
def calculate_thermal_time(T_min, T_max, T_base):
    T_min = max(T_min,0)
    T_max = max(T_max,0)
    T_opt = 26
    TD_max = 37
    T_t = 0
    for r in range(1, 9):
        f_r = (1 / 2) * (1 + cos((90 / 8) * (2 * r - 1) * pi / 180))
        T_H = max(T_min + f_r * (T_max - T_min), 0)  # Degree Celsius
        if T_H < T_opt:
            T_t += T_H - T_base
        elif T_H == T_opt:
            T_t += T_opt - T_base
        elif T_H < TD_max:
            T_t += (T_opt - T_base) * (TD_max - T_H) / (TD_max - T_opt)
        else:
            T_t += 0
    return max((1 / 8) * T_t, 0)  # Degree Celsius Days

#Rate of Change of Daylength at Emergence
def calc_RoCoDLatE(latitude, julian_day):
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
    dry_matter = pd.DataFrame({'Cohort':[0],'#Tillers':[0],'N_n':[0],'Proportion Surviving':[0],
                               'Leaf Number':[0],'Leaf Active Area':[0],'Stage':[""],'Rate of Leaf Growth':[0],'Life Span':[0]})
    dry_matter['Cohort'] = dry_matter['Cohort'].convert_dtypes(convert_integer=True)
    dry_matter['Leaf Number'] = dry_matter['Leaf Number'].convert_dtypes(convert_integer=True)

    ### Tiller and Leaf Growth Submodel ###
    #Initialising Variables
    rate_of_change_of_daylength_at_emergence = 0
    new_tillers = 0
    new_leaves = 0
    successive_leaf_thermal_time = 0
    offset_total_thermal_time = 0
    max_tillers = 0
    week_day = 0
    number_of_cohorts = 0
    #Constants for Tiller Death Proportion
    A=825
    alpha = 1.46
    beta = 2.24
    print(plant_data)
    for index,row in plant_data.iterrows():
        julian_day = row['Date'].timetuple().tm_yday
        if row['Stage'] == 'Seeding':
            continue
        elif row['Stage'] == 'Emergence':
            #Offsetting Degree Days so that we have cumulative degree days from emergence onwards to calculate new leaf growth
            if row['Stage Sum Degree Days'] == row['Daily Degree Days']:
                row['Total Degree Days'] -= (row['Total Degree Days'] - row['Daily Degree Days'])

            #Rate 
            if rate_of_change_of_daylength_at_emergence == 0:
                rate_of_change_of_daylength_at_emergence = calc_RoCoDLatE(Lat, julian_day)
                rate_of_leaf_appearance_per_degree_day = 0.025 * rate_of_change_of_daylength_at_emergence + 0.0104
                phylochron_interval = 1/rate_of_leaf_appearance_per_degree_day

            #Growing first 3 leaves
            if (row['Total Degree Days'] - successive_leaf_thermal_time) >= phylochron_interval:
                dry_matter.loc[new_leaves,'Leaf Number'] = new_leaves + 1
                dry_matter.loc[new_leaves,['Leaf Active Area','Stage','Rate of Leaf Growth']] = [0,'Grow',leaf_data[new_leaves]["max_leaf_area"] * rate_of_leaf_appearance_per_degree_day]
                successive_leaf_thermal_time = row['Total Degree Days']
                new_leaves += 1

            if dry_matter['Leaf Number'].values.max() < 3:
                continue
            else:
                week_day += 1
                new_tillers += max(row['Mean Temp'],0) * TPr * 250
                #new_tillers += calculate_thermal_time(row['Min Temp'],row['Max Temp'],T_base=1) * TPr * 250
                if week_day == 7:
                    dry_matter.loc[number_of_cohorts,['Cohort','#Tillers']] = [number_of_cohorts+1,new_tillers]
                    if number_of_cohorts == 0:
                        dry_matter.loc[number_of_cohorts,'N_n'] = 0 
                    else:
                        dry_matter.loc[number_of_cohorts,'N_n'] = dry_matter['#Tillers'][number_of_cohorts-1] + dry_matter['N_n'][number_of_cohorts-1]
                    number_of_cohorts += 1
                    week_day = 0
                    new_tillers = 0
        elif row['Stage'] == 'Double Ridge':
            #Letting last week of growth in emergence phase be counted
            if week_day > 0:
                    dry_matter.loc[number_of_cohorts,['Cohort','#Tillers']] = [number_of_cohorts+1,new_tillers]
                    if number_of_cohorts == 0:
                        dry_matter.loc[number_of_cohorts,'N_n'] = 0 
                    else:
                        dry_matter.loc[number_of_cohorts,'N_n'] = dry_matter['#Tillers'][number_of_cohorts-1] + dry_matter['N_n'][number_of_cohorts-1]
                    number_of_cohorts += 1
                    week_day = 0
                    new_tillers = 0
            for c in dry_matter['Cohort']:
                dry_matter.loc[c-1,'Proportion Surviving'] = 1 / (1 + ((row['Stage Sum Degree Days']/400) / (A/dry_matter['N_n'][c-1])**alpha)**beta)
            if (row['Total Degree Days'] - successive_leaf_thermal_time) >= phylochron_interval:
                dry_matter.loc[new_leaves,'Leaf Number'] = new_leaves + 1
                dry_matter.loc[new_leaves,['Leaf Active Area','Stage','Rate of Leaf Growth']] = [0,'Grow',leaf_data[new_leaves]["max_leaf_area"] * rate_of_leaf_appearance_per_degree_day]
                successive_leaf_thermal_time = row['Total Degree Days']
                new_leaves += 1
        if dry_matter['Rate of Leaf Growth'].values.max() != 0:
            i=0
            for area in dry_matter['Leaf Active Area']:
                if dry_matter.loc[i,'Stage'] == 'Grow':
                    if area == 0:
                        dry_matter.loc[i,'Leaf Active Area'] = max(row['Daily Degree Days'],0) * dry_matter['Rate of Leaf Growth'][i]
                    else:
                        if dry_matter.loc[i,'Leaf Active Area'] != leaf_data[i]['max_leaf_area']:
                            if dry_matter.loc[i,'Leaf Active Area'] + max(row['Daily Degree Days'],0) * dry_matter['Rate of Leaf Growth'][i] > leaf_data[i]['max_leaf_area']:
                                dry_matter.loc[i,'Leaf Active Area'] = leaf_data[i]['max_leaf_area']
                            else:
                                dry_matter.loc[i,'Leaf Active Area'] += max(row['Daily Degree Days'],0) * dry_matter['Rate of Leaf Growth'][i]
                        else:
                            dry_matter.loc[i,'Stage'] = 'Max Area'
                            if i < 9:
                                dry_matter.loc[i,'Life Span'] = 0.67 * (3.5 * phylochron_interval)
                            elif i == 9:
                                dry_matter.loc[i,'Life Span'] = 0.67 * (375)
                            else:
                                dry_matter.loc[i,'Life Span'] = 0.67 * (375 + (i-9) * 125)
                            dry_matter.loc[i,'Rate of Leaf Growth'] = 0
                elif dry_matter.loc[i,'Stage'] == 'Max Area':
                    if (dry_matter.loc[i,'Life Span'] - max(row['Daily Degree Days'],0)) > 0:
                        dry_matter.loc[i,'Life Span'] -= max(row['Daily Degree Days'],0)
                    else:
                        dry_matter.loc[i,'Stage'] = 'Decay'
                        if i < 9:
                            dry_matter.loc[i,'Life Span'] = 0.33 * (3.5 * phylochron_interval)
                        elif i == 9:
                            dry_matter.loc[i,'Life Span'] = 0.33 * (375)
                        else:
                            dry_matter.loc[i,'Life Span'] = 0.33 * (375 + (i-9) * 125)
                        dry_matter.loc[i,'Rate of Leaf Growth'] = (-1) * leaf_data[i]['max_leaf_area'] / dry_matter.loc[i,'Life Span']
                elif dry_matter.loc[i,'Stage'] == 'Decay':
                    if (dry_matter.loc[i,'Life Span'] - max(row['Daily Degree Days'],0)) > 0:
                        dry_matter.loc[i,'Life Span'] -= max(row['Daily Degree Days'],0)
                        dry_matter.loc[i,'Leaf Active Area'] += max(row['Daily Degree Days'],0) * dry_matter['Rate of Leaf Growth'][i]
                    else:
                        dry_matter.loc[i,['Leaf Active Area','Stage','Rate of Leaf Growth','Life Span']] = [0,'Dead',0,0]
                i+=1
        projected_area_factor = 1e-6  # Convert from mm^2 to m^2

        LAI_z = pd.DataFrame({'Level':[int()],'Height':[float()],'Leaf Numbers':[""],'LAI':[float()]})

        level = 1
        for leaf in dry_matter['Leaf Number'].dropna():
            print(leaf,level)
            print(leaf_data[leaf-1]['sheath_length'],'<',leaf_data[leaf]['sheath_length'],leaf_data[leaf-1]['sheath_length'] < leaf_data[leaf]['sheath_length'])
            LAI_z.loc[level,['Level','Height','LAI']] = [level,leaf_data[leaf-1]['sheath_length'],LAI_z.loc[level-1,'LAI']+dry_matter.loc[leaf-1,'Leaf Active Area']]
            if leaf_data[leaf-1]['sheath_length'] < leaf_data[leaf]['sheath_length']:
                level += 1

                

        # number_of_levels = 0
        # if dry_matter['Leaf Number'].values[-1] >= 12:
        #         number_of_levels = 5
        # else:
        #     for i in range(trunc(dry_matter['Leaf Number'].values[-1])+1):
        #         if leaf_data[i]['sheath_length'] > leaf_data[i-1]['sheath_length']:
        #             number_of_levels += 1
        # level = 0
        # i=0
        # while number_of_levels-level>=0:
        #     if dry_matter['Leaf Number'].values[-1] - i < 0:
        #         break
        #     if trunc(dry_matter['Leaf Number'].values[-1] - i) == 0:
        #         if i < 11:
        #             LAI += leaf_data[i]['leaf_area'] * (dry_matter['Leaf Number'].values[-1] - i)
        #             if leaf_data[i]['sheath_length'] > leaf_data[i-1]['sheath_length']:
        #                 level += 1
        #         else:
        #             LAI += leaf_data[11]['leaf_area'] * (dry_matter['Leaf Number'].values[-1] - i)
        #         LAI_z.loc[number_of_levels-level] = [f'Level {number_of_levels-level}',LAI]
        #         break
        #     if i > 0:
        #         if i < 11:
        #             LAI += leaf_data[i]['leaf_area']
        #             if leaf_data[i]['sheath_length'] > leaf_data[i-1]['sheath_length']:
        #                 level += 1
        #         else:
        #             LAI += leaf_data[11]['leaf_area']
        #         LAI_z.loc[number_of_levels-level] = [f'Level {number_of_levels-level}',LAI]
        #     else:
        #         LAI += leaf_data[0]['leaf_area']
        #     i+=1
        # LAI_z = LAI_z.sort_index()
        # LAI_z['LAI'] = LAI_z['LAI'].multiply(projected_area_factor)
        if len(dry_matter['Leaf Number'].values) > 1:
            print(max(dry_matter['Leaf Number'].values.dropna()))
        print(LAI_z)    
print(dry_matter)       