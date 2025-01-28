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
    dry_matter = pd.DataFrame({'Cohort':[0],'#Tillers':[0],'Total Tillers':[0],'Leaf Number':[0],'Leaf Active Area':[0],'Rate of Leaf Growth':[0]})
    dry_matter['Cohort'] = dry_matter['Cohort'].convert_dtypes(convert_integer=True)
    dry_matter['Leaf Number'] = dry_matter['Leaf Number'].convert_dtypes(convert_integer=True)

    ### Tiller and Leaf Growth Submodel ###
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
                leaf_interval = phylochron_interval

            #Growing first 3 leaves
            if (row['Total Degree Days'] - successive_leaf_thermal_time) >= phylochron_interval:
                new_leaves += 1
                dry_matter.loc[new_leaves,'Leaf Number'] = new_leaves
                dry_matter.loc[new_leaves,'Rate of Leaf Growth'] = leaf_data[new_leaves-1]["max_leaf_area"] / phylochron_interval
                successive_leaf_thermal_time = row['Total Degree Days']

            if dry_matter['Leaf Number'].values.max() < 3:
                continue
            else:
                week_day += 1
                new_tillers += max((row['Max Temp']+row['Min Temp'])/2,1) * TPr * 250
                if week_day == 7:
                    print(dry_matter['#Tillers'].values)
                    print(max(dry_matter['#Tillers'].values))
                    dry_matter.loc[number_of_cohorts,['Cohort','#Tillers','Total Tillers']] = [number_of_cohorts+1,new_tillers,max(dry_matter['Total Tillers'].values)+new_tillers if number_of_cohorts>0 else new_tillers]
                    number_of_cohorts += 1
                    week_day = 0
                    new_tillers = 0
        elif row['Stage'] == 'Double Ridge':
            if week_day > 0:
                    dry_matter.loc[number_of_cohorts,['Cohort','#Tillers','Total Tillers']] = [number_of_cohorts+1,new_tillers,max(dry_matter['Total Tillers'].values)+new_tillers if number_of_cohorts>0 else new_tillers]
                    number_of_cohorts += 1
                    week_day = 0
                    new_tillers = 0
            if (row['Total Degree Days'] - successive_leaf_thermal_time) >= phylochron_interval:
                new_leaves += 1
                dry_matter.loc[new_leaves,'Leaf Number'] = new_leaves
                dry_matter.loc[new_leaves,'Rate of Leaf Growth'] = leaf_data[new_leaves-1]["max_leaf_area"] / phylochron_interval
                successive_leaf_thermal_time = row['Total Degree Days']
            
        print(dry_matter)

print(sum(dry_matter['#Tillers']))            
    #         if max_tillers == 0:
    #             dry_matter['#Tillers'].loc[index-1] = dry_matter['#Tillers'].values[-1]
    #             max_tillers = trunc(dry_matter['#Tillers'].values[-1])+1
    #             N_n = pd.DataFrame({'N_n':[i for i in range(1,max_tillers+1)],
    #                                 'Survival Chance':[1*(dry_matter['#Tillers'].values[-1] - trunc(dry_matter['#Tillers'].values[-1])) if i == trunc(dry_matter['#Tillers'].values[-1])+1 else 1 for i in range(1,max_tillers+1)]})
    #             chance = list(int(i) for i in list('1'*max_tillers))
    #         for index2,row2 in N_n.iterrows():
    #             chance[index2] *= (1 / (1 + (((row['Stage Sum Degree Days']/400)/((A/row2['N_n'])**alpha)))**beta))
    #         print(chance)
    #         dry_matter.loc[index] = [new_tillers, new_leaves]             

    #     #N_n = pd.concat([N_n, pd.DataFrame({'Multiplier':chance})], axis=1)
    #     plant_data = pd.concat([plant_data, dry_matter], axis=1)

    #     LAI_z = pd.DataFrame({'Level':[],'LAI':[]})
    #     LAI_z['LAI'] = LAI_z['LAI'].astype(np.float64)
    #     LAI_z['Level'] = LAI_z['Level'].astype(str)
    #     LAI = float(0)

    #     number_of_levels = 0
    #     if dry_matter['#Leaves'].values[-1] >= 12:
    #             number_of_levels = 5
    #     else:
    #         for i in range(trunc(dry_matter['#Leaves'].values[-1])+1):
    #             if leaf_data[i]['sheath_length'] > leaf_data[i-1]['sheath_length']:
    #                 number_of_levels += 1
    #     level = 0
    #     i=0
    #     while number_of_levels-level>=0:
    #         if dry_matter['#Leaves'].values[-1] - i < 0:
    #             break
    #         if trunc(dry_matter['#Leaves'].values[-1] - i) == 0:
    #             if i < 11:
    #                 LAI += leaf_data[i]['max_leaf_area'] * (dry_matter['#Leaves'].values[-1] - i)
    #                 if leaf_data[i]['sheath_length'] > leaf_data[i-1]['sheath_length']:
    #                     level += 1
    #             else:
    #                 LAI += leaf_data[11]['max_leaf_area'] * (dry_matter['#Leaves'].values[-1] - i)
    #             LAI_z.loc[number_of_levels-level] = [f'Level {number_of_levels-level}',LAI]
    #             break
    #         if i > 0:
    #             if i < 11:
    #                 LAI += leaf_data[i]['max_leaf_area']
    #                 if leaf_data[i]['sheath_length'] > leaf_data[i-1]['sheath_length']:
    #                     level += 1
    #             else:
    #                 LAI += leaf_data[11]['max_leaf_area']
    #             LAI_z.loc[number_of_levels-level] = [f'Level {number_of_levels-level}',LAI]
    #         else:
    #             LAI += leaf_data[0]['max_leaf_area']
    #         i+=1
    #     LAI_z = LAI_z.sort_index()
    #     LAI_z['LAI'] = LAI_z['LAI'].multiply(projected_area_factor)
    
    #     # ### Root Growth Submodel ###
    #     # root_growth = pd.DataFrame({'Layer':[], 'Length':[], 'Weight':[]})
    #     # seminal_weight = 1.5 * (10**-4)
    #     # lateral_weight = 4 * (10**-5)
    #     # for index,row in plant_data.iterrows():
    #     #     TR = min(0.2 + 0.12 * plant_data['Mean_Temp'],0)
    #     #     if 'Seminal' in root_growth['Layer'].values(-1):
    #     #         if index > 0:
    #     #             length = TR + root_growth['Length'].values(-1)
    #     #         else: 
    #     #             length = TR
    #     #         root_growth.loc[index] = ['Seminal',length,seminal_weight*length]
    #     #     else:
    #     #         root_growth.loc[index]

    #     ### Light interception and photosynthesis submodel ###
    #     Qp_z = pd.DataFrame({'Qp_0':PAR_data.loc[julian_day-1].values[1:]})
    #     for i in range(1,number_of_levels+1):
    #         Qp_z.insert(i,f'Qp_{i}',Qp_z['Qp_0'].apply(func=(lambda x: ((x*k)/(1-m))*exp(-k*LAI_z['LAI'].values[i]))))
    # print(dry_matter)        
