from numpy import float64
import pandas as pd
import datetime
import os
from math import isnan, pi, cos, radians, sin, sqrt, exp

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\\..'))

configdir = os.path.join(project_path, 'config.txt')
with open(configdir, 'r') as config_file:
    stages = []
    for line in config_file:
        if line.__contains__(','):
            if line[:line.index(',')] == 'shoot_production_rate':
                TPr = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'Latitude':
                Lat = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'extinction_coefficient':
                k = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'leaf_transmission_coefficient':
                m = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'crop_boundary_layer_resistance':
                r_a = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'mesophyll_resistance':
                r_m = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'ambient_co2_concentration':
                C_a = float(line.strip().split(',')[1])


#Leaf data from the paper (Table 1)
leaf_data = pd.DataFrame({"Lamina Length":[125,125,125,125,125,125,245,275,310,350,395,445],
                          "Lamina Width":[6,6,6,6,6,6,11,13,15,15,15,15],
                          "Sheath Length":[0,0,0,0,0,0,120,125,125,135,145,155],
                          "Max Leaf Area":[653,653,653,653,653,979,1797,2322,3164,3481,3988,4560]})

path = os.path.join(project_path,'Data','Processed','Thermal Time')

temp_data = pd.read_csv(os.path.join(project_path, 'Data', 'Raw', 'Temperature 1978-1981.csv'))
temp_data['Date'] = pd.to_datetime(temp_data['Date'])
temp_data = temp_data.rename(columns={'Min_Temp':'Min Temp','Max_Temp':'Max Temp','Mean_Temp':'Mean Temp'})
PAR_data = pd.read_csv(os.path.join(project_path, 'Data', 'Raw', 'Average Hourly PAR.csv'),index_col=0,header=0)

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

#For testing i would recommend only have 1 file of plant data in the thermal time folder in the processed data folder
for file in os.listdir(path):
    plant_data = pd.read_csv(os.path.join(path,file))
    plant_data['Date'] = pd.to_datetime(plant_data['Date'])
    plant_data = plant_data.merge(temp_data, on='Date', how='outer')
    plant_data = plant_data.dropna().reset_index().drop(columns='index')
    
    #Setting up Dataframes
    dry_matter = pd.DataFrame({'Cohort':[int(0)],'#Tillers':[0],'N_n':[0],'Proportion Surviving':[0],
                               'Leaf Number':[int(0)],'Leaf Active Area':[0],'Stage':[""],'Rate of Leaf Growth':[float(0)],'Life Span':[0]})
    dry_matter['Cohort'] = dry_matter['Cohort'].convert_dtypes(convert_integer=True)
    dry_matter['Leaf Number'] = dry_matter['Leaf Number'].convert_dtypes(convert_integer=True)

    LAI_z = pd.DataFrame({'Level':[0],'Height':[0],'LAI':[0]})
    LAI_z['Level'] = LAI_z['Level'].convert_dtypes(convert_integer=True)
    LAI_z['LAI'] = LAI_z['LAI'].astype(dtype=float64)

    Photosynthesis = pd.DataFrame({'Level (z)':[],'Qp':[]})
    Photosynthesis['Level (z)'] = Photosynthesis['Level (z)'].convert_dtypes(convert_integer=True)
    Photosynthesis['Qp'] = Photosynthesis['Qp'].astype(dtype=float64)

    hourly_temp = pd.DataFrame()

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

    for index,row in plant_data.iterrows():
        julian_day = row['Date'].timetuple().tm_yday
        if row['Stage'] == 'Seeding':
            continue
        elif row['Stage'] == 'Emergence':
            #Offsetting Degree Days so that we have cumulative degree days from emergence onwards to calculate new leaf growth
            if row['Stage Sum Degree Days'] == row['Daily Degree Days']:
                row['Total Degree Days'] -= (row['Total Degree Days'] - row['Daily Degree Days'])

            #Rate of new leaves
            if rate_of_change_of_daylength_at_emergence == 0:
                rate_of_change_of_daylength_at_emergence = calc_RoCoDLatE(Lat, julian_day)
                rate_of_leaf_appearance_per_degree_day = 0.025 * rate_of_change_of_daylength_at_emergence + 0.0104
                phylochron_interval = 1/rate_of_leaf_appearance_per_degree_day

            #Growing first 3 leaves
            if (row['Total Degree Days'] - successive_leaf_thermal_time) >= phylochron_interval:
                dry_matter.loc[new_leaves,'Leaf Number'] = new_leaves + 1
                dry_matter.loc[new_leaves,['Leaf Active Area','Stage','Rate of Leaf Growth']] = [0,'Grow',leaf_data.loc[new_leaves,"Max Leaf Area"] * rate_of_leaf_appearance_per_degree_day]
                successive_leaf_thermal_time = row['Total Degree Days']
                new_leaves += 1

            if dry_matter['Leaf Number'].values.max() < 3:
                continue
            else:
                week_day += 1
                #new_tillers += max(row['Mean Temp'],0) * TPr * 250
                new_tillers += calculate_thermal_time(row['Min Temp'],row['Max Temp'],T_base=1) * TPr * 250
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
            #Calculating the proportion surviving in each Cohort
            for c in dry_matter['Cohort'].dropna():
                if c == 1:
                    dry_matter.loc[c-1,'Proportion Surviving'] = 1
                else:
                    dry_matter.loc[c-1,'Proportion Surviving'] = 1 / (1 + ((row['Stage Sum Degree Days']/400) / (A/dry_matter['N_n'][c-1])**alpha)**beta)
            #Grow Leaves
            if (row['Total Degree Days'] - successive_leaf_thermal_time) >= phylochron_interval:
                dry_matter.loc[new_leaves,'Leaf Number'] = new_leaves + 1
                dry_matter.loc[new_leaves,['Leaf Active Area','Stage','Rate of Leaf Growth']] = [0,'Grow',leaf_data.loc[new_leaves,"Max Leaf Area"] * rate_of_leaf_appearance_per_degree_day]
                successive_leaf_thermal_time = row['Total Degree Days']
                new_leaves += 1
        #Leaf Growth
        if dry_matter['Rate of Leaf Growth'].values.max() != 0:
            i=0
            for area in dry_matter['Leaf Active Area']:
                #Growth Stage
                if dry_matter.loc[i,'Stage'] == 'Grow':
                    if area == 0:
                        dry_matter.loc[i,'Leaf Active Area'] = max(row['Daily Degree Days'],0) * dry_matter['Rate of Leaf Growth'][i]
                    else:
                        if dry_matter.loc[i,'Leaf Active Area'] != leaf_data.loc[i,'Max Leaf Area']:
                            if dry_matter.loc[i,'Leaf Active Area'] + max(row['Daily Degree Days'],0) * dry_matter['Rate of Leaf Growth'][i] > leaf_data.loc[i,'Max Leaf Area']:
                                dry_matter.loc[i,'Leaf Active Area'] = leaf_data.loc[i,'Max Leaf Area']
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
                #Max Area Stage
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
                        dry_matter.loc[i,'Rate of Leaf Growth'] = (-1) * leaf_data.loc[i,'Max Leaf Area'] / dry_matter.loc[i,'Life Span']
                #Decay Stage
                elif dry_matter.loc[i,'Stage'] == 'Decay':
                    if (dry_matter.loc[i,'Life Span'] - max(row['Daily Degree Days'],0)) > 0:
                        dry_matter.loc[i,'Life Span'] -= max(row['Daily Degree Days'],0)
                        dry_matter.loc[i,'Leaf Active Area'] += max(row['Daily Degree Days'],0) * dry_matter['Rate of Leaf Growth'][i]
                    else:
                        dry_matter.loc[i,['Leaf Active Area','Stage','Rate of Leaf Growth','Life Span']] = [0,'Dead',0,0]
                i+=1
        #Leaf Area Index

        level = 1
        for leaf in dry_matter['Leaf Number'].values.dropna():
            LAI_z.loc[level,'Level'] = level
            LAI_z.loc[level,'Height'] = round(leaf_data.loc[leaf-1,'Sheath Length'] * 1e-3 * 250, 8)
            LAI_z.loc[level,'LAI'] = round(LAI_z.loc[level-1,'LAI']+dry_matter.loc[leaf-1,'Leaf Active Area'] * 1e-6 * 250, 8)

            if leaf_data.loc[leaf-1,'Sheath Length'] < leaf_data.loc[leaf,'Sheath Length']:
                level += 1
            elif level == 0:
                level += 1
                
        #         ### Root Growth Submodel ###
        # root_growth = pd.DataFrame({'Layer':[], 'Length':[], 'Weight':[]})
        # seminal_weight = 1.5 * (10**-4)
        # lateral_weight = 4 * (10**-5)
        # TR = min(0.2 + 0.12 * row['Mean Temp'],0)
        # if 'Seminal' in root_growth['Layer'].values[-1]:
        #     if index > 0:
        #         length = TR + root_growth['Length'].values(-1)
        #     else: 
        #         length = TR
        #     root_growth.loc[index] = ['Seminal',length,seminal_weight*length]
        # else:
        #     root_growth.loc[index]

        ## Light interception and photosynthesis submodel ###
        for i in LAI_z['Level']:
            Photosynthesis.loc[i,'Level (z)'] = i
            H = 0
            R = 0
            hour = 0
            sunrise = 0
            for j in PAR_data.loc[julian_day]:
                if j!= 0 and sunrise == 0:
                    sunrise = hour
                elif j==0 and sunrise != 0:
                    sunset = hour
                else:
                    hour += 1
            hour = 0
            for j in PAR_data.loc[julian_day]:
                #Counting number of daylight hours
                if j != 0:
                    H += 1
                #Qp is the intensity of PAR at a given layer
                try:
                    if LAI_z.loc[i,'Level'] != 0:
                        Qp = ((Qp*k)/(1-m))*exp((-k)*(LAI_z.loc[i-1,'LAI']))
                    else:
                        Qp = j
                except NameError:                 
                    if LAI_z.loc[i,'Level'] == 0:
                        Qp = j
                if Qp != 0:
                    #r_s is the stomatal resistance
                    r_s = 1.56 * 75*(1+(100/Qp)) #*(1-0.3*D) is usually included
                    #where D is the vapour pressure deficity however this crop is considered to be free from water stress
                    #r_p is the total physical resistance
                    r_p = r_a + r_s + r_m
                    #P_max is hte maximum photosynthesis rate
                    P_max = 0.995*(C_a/r_p)
                    #I do not have hourly temperature but I will assume it follows a trigonometric pattern between peaks
                    #with max temp at 13:00 and min temp in hour last hour without sunlight
                else:
                    P_max = 0
                if hour < sunrise:
                    T_h = ((plant_data.loc[julian_day-1,'Max Temp']-row['Min Temp'])/2)*cos(((hour-12)*pi)/(24-(13-sunrise))) + ((plant_data.loc[julian_day-1,'Max Temp']+row['Min Temp'])/2)
                elif sunrise <= hour and hour <= 13:
                    T_h = ((row['Max Temp']-row['Min Temp'])/2)*sin(((hour-0.5)*pi)/(13-sunrise)) + ((row['Max Temp']+row['Min Temp'])/2)
                else:
                    T_h = ((row['Max Temp']-plant_data.loc[julian_day+1,'Min Temp'])/2)*cos(((hour-12)*pi)/(24-(13-sunrise))) + ((row['Max Temp']+plant_data.loc[julian_day+1,'Min Temp'])/2)
                #print(f"{hour}:{T_h}")
                hourly_temp.loc[julian_day,hour] = T_h
                # #P_m is the temperature-dependent maximum photosynthetic rate
                # P_m = (0.044*6*(10**9)*)
                # #R is the Total respiration in each day
                # R += 0.65*grc*

                hour += 1
        print(Photosynthesis)
        print(dry_matter)
        print(hourly_temp)
        
        ###Testing number of tillers grown###
        if row['Stage'] == 'Double Ridge':
            raise KeyboardInterrupt

    print(file)
    print(LAI_z)    
    print(dry_matter)       