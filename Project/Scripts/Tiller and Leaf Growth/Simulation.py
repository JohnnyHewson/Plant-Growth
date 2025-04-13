from numpy import float64
import pandas as pd
import datetime
import os
from math import acos, pi, cos, radians, sin, sqrt, exp, tan
import matplotlib.pyplot as plt

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\\..'))

configdir = os.path.join(project_path, 'config.txt')
with open(configdir, 'r') as config_file:
    stages = []
    for line in config_file:
        if line.__contains__(','):
            if line[:line.index(',')] in ['Seeding',
                                          'Emergence',
                                          'Double Ridge',
                                          'Anthesis',
                                          'Maturity']:
                stage_name, T_base, degree_days = line.strip().split(',')
                stages.append({
                    'Stage': stage_name,
                    'T_Base': float(T_base),
                    'Degree_Days': float(degree_days)
                    })
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
            elif line[:line.index(',')] == 'activation_energy_of_the_electron_transport_system':
                dH_1 = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'denaturation_energy_of_the_electron_transport_system':
                dH_2 = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'entropy_change_on_denaturation_of_the_electron_transport_system':
                dS = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'growth_respiration_coefficient':
                grc = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'maintenance_respiration_coefficient_emerge_to_anthesis':
                mrc_e2a = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'maintenance_respiration_coefficient_anthesis_to_maturity':
                mrc_a2m = float(line.strip().split(',')[1])
            elif line[:line.index(',')] == 'emergence_to_double_ridge':
                assimDistE2DR = tuple(line.strip().split(',')[1:5])
            elif line[:line.index(',')] == 'double_ridge_to_beginning_ear_growth':
                assimDistDR2BEG = tuple(line.strip().split(',')[1:5])
            elif line[:line.index(',')] == 'beginning_ear_growth_to_anthesis':
                assimDistBEG2A = tuple(line.strip().split(',')[1:5])


#Leaf data from the paper (Table 1)
leaf_data = pd.DataFrame({"Lamina Length":[125,125,125,125,125,125,245,275,310,350,395,445],
                          "Sheath Length":[0,0,0,0,0,0,120,125,125,135,145,155],
                          "Max Leaf Area":[653,653,653,653,653,979,1797,2322,3164,3481,3988,4560],
                          "Layer":[1,1,1,1,1,1,2,3,3,4,5,6]})

path = os.path.join(project_path,'Data','Processed','Thermal Time')

temp_data = pd.read_csv(os.path.join(project_path, 'Data', 'Raw', 'Temperature 1978-1981.csv'))
temp_data['Date'] = pd.to_datetime(temp_data['Date'])
temp_data = temp_data.rename(columns={'Min_Temp':'Min Temp','Max_Temp':'Max Temp','Mean_Temp':'Mean Temp'})
PAR_data = pd.read_csv(os.path.join(project_path, 'Data', 'Raw', 'Average Hourly PAR.csv'),index_col=0,header=0)

#Overall Dataframe
Results = pd.DataFrame(index = ['Top Weight, Anthesis (g/m^2)',
                                'Top Weight, Maturity (g/m^2)',
                                'Grain Yield (g/m^2)',
                                'Harvest Index (%)',
                                'No. of Grain per m^2 (10^3)',
                                'No. of Ears per m^2',
                                'No. Grains Per Ear',
                                'Grain Weight (mg/grain)',
                                'Grain Pool, Anthesis (g/m^2)',
                                'Grain Pool, Maturity (g/m^2)',
                                'Root Weight, Anthesis (g/m^2)'])

#Calculates T_H for Thermal Time and Vernalized Degree Days Calculations
def calc_T_H(T_max, T_min, r):
    f_r = (1 / 2) * (1 + cos((90 / 8) * (2 * r - 1) * pi / 180))
    T_H = max(T_min + f_r * (T_max - T_min), 0)  # Degree Celsius

    return T_H

#Calculates Un/affected Thermal Time (Celcius)
def calc_thermal_time(T_min, T_max, VDD, declination, latitude, current_stage, stages, affected):
    T_min = max(T_min,0)
    T_max = max(T_max,0)
    T_opt = 26
    TD_max = 37
    if affected:
        T_base = next(dic['T_Base'] for dic in stages if current_stage == dic['Stage'])
    else:
        T_base = 0

    T_t = 0
    for r in range(1, 9):
        T_H = calc_T_H(T_max, T_min, r)
        if T_H < T_opt:
            T_t += T_H - T_base
        elif T_H == T_opt:
            T_t += T_opt - T_base
        elif T_opt < TD_max:
            T_t += (T_opt - T_base) * (TD_max - T_H) / (TD_max - T_opt)
        else:
            T_t += 0
    T_t = max((1 / 8) * T_t, 0)

    if affected:
        return T_t * calc_FP(declination, latitude, current_stage) * calc_FV(current_stage, VDD)
    else:
        return T_t


def calc_FP(declination, latitude, current_stage):
    #FP only affects thermal time during Emergence and Double Ridge
    if current_stage not in ["Emergence", "Double Ridge"]:
        return 1
    elif current_stage == "Double Ridge":
        P_base = 7
    elif current_stage == "Emergence":
        P_base = 0
    P_opt = 20
    latitude = radians(latitude)

    D = (-0.10453) / (cos(latitude) * cos(declination))
    P_R = acos(D - (tan(latitude) * tan(declination)))
    P_H = P_R * (24/pi)
    FP = (P_H - P_base) / (P_opt - P_base)

    return min(max(FP,0),1)

def calc_FV(current_stage, VDD):
    if current_stage != "Emergence":
        return 1
    else:
        V_base = 8
        V_sat = 33

        FV = (VDD - V_base) / (V_sat - V_base)

        return min(max(FV,0),1)

def calc_VDD(T_max, T_min, VDD):
    T_min = max(T_min,0)
    T_max = max(T_max,0)

    if T_max > 30:
        VDD /= 2

    V_eff = 0
    for r in range(1,9):
        T_H = calc_T_H(T_max, T_min, r)
        if -4 <= T_H and T_H < 3:
            V_eff += (T_H + 4) / 7
        elif 3 <= T_H and T_H < 10:
            V_eff += 1
        elif 10 <= T_H and T_H < 17:
            V_eff += (17 - T_H) / 7

    VDD += (1/8) * V_eff

    return VDD

#Calculate Declination Angle
def calc_dec(julian_day):
    tilt_of_earth = 23.44 * (pi / 180) #Converts to radians
    radians_per_day = 2*pi * (1 / 365.25)
    declination = -tilt_of_earth * cos(radians_per_day * (julian_day + 10))
    
    return declination

#Rate of Change of Daylength at Emergence
def calc_RoCoDLatE(latitude, declination):
    tilt_of_earth = 23.44 * (pi / 180) #Converts to radians
    radians_per_day = 2*pi * (1 / 365.25)

    latitude = radians(latitude)
    declination_rate = (tilt_of_earth * radians_per_day * sin(radians_per_day * (julian_day + 10)))

    numerator = cos(declination) * sin(latitude)
    denominator = sqrt((1 - (sin(latitude) * sin(declination)))**2)
    
    if denominator == 0:
        return 0.0
    else:
        return (24 / pi) * (numerator / denominator) * declination_rate

def plot_tillers(values, labels):
    plt.figure(figsize=(8, 5))
    
    for i in range(len(labels)):  # 9 running lines
        x_values = [x_value for x_value, _ in values]
        y_values = [y_series.iloc[i] for _, y_series in values]
        label_name = labels.iloc[i] if i < len(labels) else f'Position {i}'
        plt.plot(x_values, y_values, marker='', linestyle='-', linewidth=2, label=label_name)
    
    plt.xlabel('Thermal Time')
    plt.ylabel('Proportion of Cohort Surviving')
    plt.title('Proportion of Cohort Surviving vs. Thermal Time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_peakLAI(LAIGraph):
    plt.plot([x for x,_ in LAIGraph],[y for _,y in LAIGraph])
    plt.xlabel('Thermal Time')
    plt.ylabel('Peak LAI')
    plt.title('Peak LAI vs. Thermal Time')
    plt.legend()
    plt.show()

#For testing i would recommend only have 1 file of plant data in the thermal time folder in the processed data folder
seeding_dates = pd.read_csv(os.path.join(project_path, 'Data', 'Raw','Seedings Dates.csv'),index_col=0,header=0)
for Plant_ID, date in seeding_dates.iterrows():
    plant_data = pd.DataFrame()
    plant_data['Date'] = pd.date_range(start=date["Seeding Date"], end=datetime.datetime((datetime.datetime.strptime(date["Seeding Date"], "%Y-%m-%d").year + 1),8,31)) #Set date range from day of seeding to end of growing year
    plant_data = plant_data.merge(temp_data[(temp_data['Date'] >= plant_data['Date'].iloc[0]) &
                                            (temp_data['Date'] <= plant_data['Date'].iloc[-1])], on='Date', how='outer')
    plant_data = plant_data.dropna().reset_index().drop(columns='index')

    Results[Plant_ID] = [0,0,0,0,0,0,0,0,0,0,0]
    #Setting up Dataframes
    dry_matter = pd.DataFrame({'Cohort':[int(0)],'#Tillers':[0],'N_n':[0],'Proportion Surviving':[0],
                               'Leaf Number':[int(0)],'Leaf Active Area':[float64(0)],'Stage':[""],'Rate of Leaf Growth':[float(0)],'Life Span':[0]})
    dry_matter['Cohort'] = dry_matter['Cohort'].convert_dtypes(convert_integer=True)
    dry_matter['Leaf Number'] = dry_matter['Leaf Number'].convert_dtypes(convert_integer=True)

    peak_tillers = pd.Series()

    LAI_z = pd.DataFrame({'Level':[int(0)],'Height':[float64(0)],'LAI':[float64(0)]})
    LAI_z['Level'] = LAI_z['Level'].convert_dtypes(convert_integer=True)

    Photosynthesis = pd.DataFrame({'Level (z)':[],'Qp':[]})
    Photosynthesis['Level (z)'] = Photosynthesis['Level (z)'].convert_dtypes(convert_integer=True)
    Photosynthesis['Qp'] = Photosynthesis['Qp'].astype(dtype=float64)

    weightDistribution = pd.DataFrame([[0,0,0,0,0,0,0]],
                                      columns=["Growth Respiration","Maintenance Respiration","Grain","Ears","Photosynthate Pool","Stems and Leaves","Roots"])
    weightDistribution.index.name = "Jday"

    Assimilate_Stage_Distribution = pd.DataFrame([(0,0,0,0),assimDistE2DR,assimDistDR2BEG,assimDistBEG2A,(0,0,0,1)],columns=['Root','Leaves','Stems','Ears'],
                                                 index=['No Assimilate','Emergence','DR Pre Grain','DR Grain','Anthesis'],dtype=float64)

    ### Tiller and Leaf Growth Submodel ###
    #Initialising Variables
    stage_number = 0
    stage_Tt = 0
    Total_Tt = 0
    VDD = 0
    rate_of_change_of_daylength_at_emergence = 0
    phylochron_interval = 0
    new_tillers = 0
    new_leaves = 0
    successive_leaf_thermal_time = 0
    offset_total_thermal_time = 0
    max_tillers = 0
    week_day = 0
    number_of_cohorts = 0
    earGrowth = False
    assimilatePool = 0
    NnPeak = pd.Series()
    peakLAI = 0
    #Constants for Tiller Death Proportion
    A=825
    alpha = 1.46
    beta = 2.24

    tillerData = []
    tillerSurvival = []
    LAIGraph = []
    
    for index,row in plant_data.iterrows():
        julian_day = row['Date'].timetuple().tm_yday
        #Assign Stage
        row['Stage'] = stages[stage_number]['Stage']
        plant_data.loc[index, 'Stage'] = row['Stage']
        ###Calculating Degree Days per day, per stage and overall
        #Per Day
        row['Daily Degree Days'] = calc_thermal_time(row['Max Temp'], row['Min Temp'], VDD, calc_dec(julian_day), Lat, row['Stage'], stages, affected=True)
        row['Daily Unaffected Thermal Time'] = calc_thermal_time(row['Max Temp'], row['Min Temp'], VDD, calc_dec(julian_day), Lat, row['Stage'], stages, affected=False)
        if 'Daily Degree Days' in plant_data:
            if index > 1 and row['Stage'] != plant_data.loc[index - 1, 'Stage']:
                #Add the leftover degree days from the last stage if the stage changes
                row['Daily Degree Days'] += (plant_data.loc[index - 2, 'Stage Sum Degree Days'] + plant_data.loc[index - 1, 'Daily Degree Days'] - stages[stage_number - 1]['Degree_Days'])
        plant_data.loc[index, 'Daily Degree Days'] = row['Daily Degree Days']
        plant_data.loc[index, 'Daily Unaffected Thermal Time'] = row['Daily Unaffected Thermal Time']

        #Per Stage
        if stage_Tt + row['Daily Degree Days'] >= stages[stage_number]['Degree_Days'] and row['Stage'] != stages[-1]['Stage']: #Maturity is the last stage
            plant_data.loc[index, 'Stage Sum Degree Days'] = stages[stage_number]['Degree_Days']
            stage_number += 1
            stage_Tt = 0
        else:
            stage_Tt += row['Daily Degree Days']
            plant_data.loc[index, 'Stage Sum Degree Days'] = stage_Tt
        row['Stage Sum Degree Days'] = plant_data.loc[index, 'Stage Sum Degree Days']

        # Step 4: Calculate Total Degree Days
        Total_Tt += row['Daily Degree Days']
        plant_data.loc[index, 'Total Degree Days'] = Total_Tt
        if index > 0:
            if plant_data.loc[index, 'Total Degree Days'] < plant_data.loc[index-1, 'Total Degree Days']:
                raise #302
        row['Total Degree Days'] = plant_data.loc[index, 'Total Degree Days']

        #Calculating VDD as its cumulative starting from germination
        VDD = calc_VDD(row['Max Temp'], row['Min Temp'], VDD)

        #Starting stage based effects
        if row['Stage'] == 'Seeding':
            offset_total_thermal_time += row['Daily Degree Days']
            continue
        elif row['Stage'] == 'Emergence':
            offset_total_thermal_time += row['Daily Degree Days']
            #Rate of new leaves
            if rate_of_change_of_daylength_at_emergence == 0:
                rate_of_change_of_daylength_at_emergence = calc_RoCoDLatE(Lat, calc_dec(julian_day))
                rate_of_leaf_appearance_per_degree_day = (0.025 * rate_of_change_of_daylength_at_emergence) + 0.0104
                phylochron_interval = 1/rate_of_leaf_appearance_per_degree_day

            #Growing leaves
            if (offset_total_thermal_time - successive_leaf_thermal_time) >= phylochron_interval:
                dry_matter.loc[new_leaves,'Leaf Number'] = new_leaves + 1
                dry_matter.loc[new_leaves,['Leaf Active Area','Stage','Rate of Leaf Growth']] = [0,'Grow',leaf_data.loc[new_leaves if new_leaves < 11 else 11,"Max Leaf Area"] * rate_of_leaf_appearance_per_degree_day/1.8] #Leaves above 12 have the same dimensions as leaf 12 (which has index 11)
                successive_leaf_thermal_time = offset_total_thermal_time
                new_leaves += 1
            
            #Growing tillers once there are 3 leaves
            if dry_matter['Leaf Number'].values.max() >= 3:
                week_day += 1
                new_tillers += max(row['Mean Temp'],0) * TPr * 250
                print(new_tillers)
                dry_matter.loc[number_of_cohorts,['Cohort','#Tillers']] = [number_of_cohorts+1,new_tillers+dry_matter.loc[number_of_cohorts-1,'N_n'] if number_of_cohorts > 0 else new_tillers]
                if number_of_cohorts == 0:
                    dry_matter.loc[number_of_cohorts,'N_n'] = 0
                else:
                    dry_matter.loc[number_of_cohorts,'N_n'] = dry_matter.loc[number_of_cohorts-1,'#Tillers'] + dry_matter.loc[number_of_cohorts-1,'N_n']
                if week_day == 7:
                    number_of_cohorts += 1
                    week_day = 0
                    #new_tillers = 0
        elif row['Stage'] == 'Double Ridge':
            offset_total_thermal_time += row['Daily Degree Days']
            #Calculating the proportion surviving in each Cohort

            # #adding 1600 tiller cohort for testing
            # try:
            #     dry_matter.loc[number_of_cohorts,['Cohort','#Tillers']]
            # except:
            #     dry_matter.loc[number_of_cohorts,['Cohort','#Tillers']] = [number_of_cohorts+1, 1600-dry_matter.loc[number_of_cohorts-1,'N_n']]
            #     dry_matter.loc[number_of_cohorts,'N_n'] = 1600
            
            if peak_tillers.empty:
                peak_tillers = dry_matter['#Tillers']

            for c in dry_matter['Cohort'].dropna():
                if c == 1:
                    dry_matter.loc[c-1,'Proportion Surviving'] = 1
                elif c > 1:
                    dry_matter.loc[c-1,'Proportion Surviving'] = 1 / (1 + ((row['Stage Sum Degree Days']/400) / (A/dry_matter.loc[c-1,'N_n'])**alpha)**beta)
            
            dry_matter['#Tillers'] = peak_tillers.mul(dry_matter['Proportion Surviving'])
            #Grow Leaves
            if (offset_total_thermal_time - successive_leaf_thermal_time) >= phylochron_interval:
                dry_matter.loc[new_leaves,'Leaf Number'] = new_leaves + 1
                dry_matter.loc[new_leaves,['Leaf Active Area','Stage','Rate of Leaf Growth']] = [0,'Grow',leaf_data.loc[new_leaves if new_leaves < 11 else 11,"Max Leaf Area"] * rate_of_leaf_appearance_per_degree_day/1.8] #Leaves above 12 have the same dimensions as leaf 12 (which has index 11)
                successive_leaf_thermal_time = offset_total_thermal_time
                new_leaves += 1
        #Leaf Growth
        if dry_matter['Rate of Leaf Growth'].values.max() != 0:
            i=0
            for area in dry_matter['Leaf Active Area']:
                #Growth Stage
                if dry_matter.loc[i,'Stage'] == 'Grow':
                    if area == 0:
                        dry_matter.loc[i,'Leaf Active Area'] = max(row['Daily Degree Days'],0) * dry_matter.loc[i,'Rate of Leaf Growth']
                    else:
                        if dry_matter.loc[i,'Leaf Active Area'] != leaf_data.loc[i if i < 11 else 11,'Max Leaf Area']:
                            if dry_matter.loc[i,'Leaf Active Area'] + max(row['Daily Degree Days'],0) * dry_matter.loc[i,'Rate of Leaf Growth'] > leaf_data.loc[i if i < 11 else 11,'Max Leaf Area']:
                                dry_matter.loc[i,'Leaf Active Area'] = leaf_data.loc[i if i < 11 else 11,'Max Leaf Area']
                            else:
                                dry_matter.loc[i,'Leaf Active Area'] += max(row['Daily Degree Days'],0) * dry_matter.loc[i,'Rate of Leaf Growth']
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
                        dry_matter.loc[i,'Leaf Active Area'] += max(row['Daily Degree Days'],0) * dry_matter.loc[i,'Rate of Leaf Growth']
                    else:
                        dry_matter.loc[i,['Leaf Active Area','Stage','Rate of Leaf Growth','Life Span']] = [0,'Dead',0,0]
                i+=1
        #Leaf Area Index
        if dry_matter['Leaf Number'].max() > 12: #Leaves 13 and above have identical dimensions to leaf 12
            leaf_data = pd.concat([leaf_data, pd.DataFrame([leaf_data.iloc[11]] * (dry_matter['Leaf Number'].max() - 12))], ignore_index=True)

        if pd.notna(dry_matter.loc[0,'Leaf Number']):
            #Get the layers that the leaves fulfill
            unique_layers = leaf_data.loc[:dry_matter['Leaf Number'].last_valid_index(), 'Layer'].unique()
            #Create a dataframe of the new info
            new_rows = pd.DataFrame({'Level': unique_layers})
            new_rows['Height'] = leaf_data['Sheath Length'].unique()[:len(unique_layers)] * 1e-3 #Convert Sheath Length from mm to m
            #This line creates a series of the cummulative sum of the current leaf areas, then using the index that it gets from the last row of the leaf data dataframe for each layer, it gets the cummulative leaf area up to an including each level
            new_rows['LAI'] = dry_matter['Leaf Active Area'].cumsum()[(min(dry_matter['Leaf Number'].last_valid_index(),i) for i in leaf_data[:dry_matter['Leaf Number'].max()].groupby('Layer').tail(1).index)].reset_index(drop=True) * 1e-6 * 250 #Convert Leaf Area from mm^2 to m^2 and 250 because paper plants/m^2
            #Re-add canopy layer
            new_rows.index += 1
            new_rows.loc[0] = [0,0,0]
            #Overwrite LAI_z
            LAI_z = new_rows.sort_index()

        ## Light interception and photosynthesis submodel ###
        P_g_h = pd.Series(0,range(0,24),dtype=float64)
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
                    if i != 0:
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
                    #P_max is the maximum photosynthesis rate
                    P_max = 0.995*(C_a/r_p)
                    #I do not have hourly temperature but I will assume it follows a trigonometric pattern between peaks
                    #with max temp at 13:00 and min temp in hour last hour without sunlight
                else:
                    P_max = 0
                if hour < sunrise:
                    T_h = ((plant_data.loc[index-1,'Max Temp']-row['Min Temp'])/2)*cos(((hour-12)*pi)/(24-(13-sunrise))) + ((plant_data.loc[index-1,'Max Temp']+row['Min Temp'])/2)
                elif sunrise <= hour and hour <= 13:
                    T_h = ((row['Max Temp']-row['Min Temp'])/2)*sin(((hour-0.5)*pi)/(13-sunrise)) + ((row['Max Temp']+row['Min Temp'])/2)
                else:
                    T_h = ((row['Max Temp']-plant_data.loc[index+1 if index+1 < len(plant_data['Date']) else index,'Min Temp'])/2)*cos(((hour-12)*pi)/(24-(13-sunrise))) + ((row['Max Temp']+plant_data.loc[index+1 if index+1 < len(plant_data['Date']) else index,'Min Temp'])/2)
                #P_g is summed of daylight hours
                if Qp != 0:
                    #P_m is the temperature-dependent maximum photosynthetic rate
                    P_m = (0.044*6*(10**9)*(T_h+273.15)*exp((-1*dH_1)/(1.987*(T_h+273.15))))/(1+(exp((-1*dH_2)/(1.987*(T_h+273.15))))*(exp(dS/1.987)))
                    #P_g is the rate of photosynthesis, this equation uses the temperature corrected quadratic equation
                    P_g_a=(0.995/P_max)*((1/alpha*Qp)+(1/P_m))
                    P_g_b=(-1*((1/alpha*Qp)+(1/P_m)+(1/P_max)))
                    P_g_c=1
                    #Taking positive P_g ###assumption
                    P_g = max((((-1*P_g_b)+(sqrt((P_g_b**2)-(4*P_g_a*P_g_c))))/(2*P_g_a)),(((-1*P_g_b)-(sqrt((P_g_b**2)-(4*P_g_a*P_g_c))))/(2*P_g_a)))
                else:
                    P_g = 0
                #Summing P_g for each layer for daily gross total
                P_g_h[hour]+=P_g*3600/1000 #Converting from per mg of CO2/sec to g of CO2/hour
                hour += 1
        if row['Stage'] in ['Emergence','Double Ridge']:
            mrc = mrc_e2a
        elif row['Stage'] in ['Anthesis','Maturity']:
            mrc = mrc_a2m
        weight = sum(sum(weightDistribution.fillna(0).values))
        growthRespiration = 0.65*grc*sum(P_g_h)
        maintenanceRespiration = weight*mrc*(2**(0.05*(row['Max Temp']+row['Min Temp'])))
        netRespiration = growthRespiration + maintenanceRespiration
        weightDistribution.loc[index,['Maintenance Respiration','Growth Respiration']] = [(0.65*grc*sum(P_g_h)),weight*mrc*(2**(0.05*(row['Max Temp']+row['Min Temp'])))]
        
                        ###Dry-matter partitioning and grain growth submodel###
        netAssimilate = (0.65 * sum(P_g_h)) - netRespiration
        
        #Setting and calculating the assimilate distribution as well as some variables that are measured at the start of the stage for the results dataframe
        if row['Stage'] == "Seeding":
            Assimilate_Distribution = Assimilate_Stage_Distribution.loc['No Assimilate']
        elif row['Stage'] == "Emergence":
            Assimilate_Distribution = Assimilate_Stage_Distribution.loc['Emergence']
        elif row['Stage'] == "Double Ridge":
            if row['Stage Sum Degree Days'] >= 200:
                earGrowth = True
                Assimilate_Distribution = Assimilate_Stage_Distribution.loc['DR Grain']
                weightDistribution.loc[index,'Ears'] = Assimilate_Distribution['Ears'] * netAssimilate
            else:
                Assimilate_Distribution = Assimilate_Stage_Distribution.loc['DR Pre Grain']
        elif row['Stage'] == "Anthesis":
            if Results.fillna(0).loc['No. Grains Per Ear',Plant_ID] == 0:
                Results.loc['No. Grains Per Ear',Plant_ID] = ((weightDistribution.loc[index,'Ears'] / 10**-2) / (dry_matter.loc[dry_matter['#Tillers'].last_valid_index(),'#Tillers'])) #10**-2 to convert 10mg to g
            Assimilate_Distribution = Assimilate_Stage_Distribution.loc['Anthesis']
            if assimilatePool == 0:
                assimilatePool = 0.3 * sum(weightDistribution["Stems and Leaves"].fillna(0).values)
            if row['Stage Sum Degree Days'] <= 55:
                assimilatePool += netAssimilate
            elif row['Stage Sum Degree Days'] <= 295:
                G_Max = ((0.045 * (row['Max Temp'] + row['Min Temp'])) / 2) + 0.4
                weightDistribution.loc[index,'Grain'] = G_Max if netAssimilate + assimilatePool > G_Max else netAssimilate + assimilatePool
            elif row['Stage Sum Degree Days'] <= 350:
                if netAssimilate < 0:
                    weightDistribution.loc[index,'Grain'] = netAssimilate

                        ### Root Growth Submodel ###
        #This submodel acts as a assimilate sink so i chose not to properly code it as it is pointless to do so
        root_assimilate = netAssimilate * Assimilate_Distribution['Root']
        TR = max(0.2 + 0.12 * row['Mean Temp'],0)
        seminalSpecificWeight = 4 * (10**(-5))
        lateralSpecificWeight = 1.5 * (10**(-4))
        seminalAssimilate = TR * (5 * seminalSpecificWeight) #5 seminal roots
        lateralAssimilate = max(root_assimilate - seminalAssimilate,0)
        seminalExtension = TR*seminalAssimilate
        lateralExtension = lateralAssimilate*lateralSpecificWeight
        seminalWeight = seminalSpecificWeight*seminalExtension
        lateralWeight = lateralSpecificWeight*lateralExtension
        root_weight = seminalWeight + lateralWeight
        weightDistribution.loc[index,"Root"] = root_weight

        #Results Dataframe Collection
        if row['Stage'] == 'Anthesis':
            if weight > Results.loc['Top Weight, Anthesis (g/m^2)',Plant_ID]:
                Results.loc['Top Weight, Anthesis (g/m^2)',Plant_ID] = weight
            if Results.loc['Grain Pool, Anthesis (g/m^2)',Plant_ID] == 0:
                Results.loc['Grain Pool, Anthesis (g/m^2)',Plant_ID] = assimilatePool
                Results.loc['Root Weight, Anthesis (g/m^2)'] = weightDistribution['Roots']
        elif row['Stage'] == 'Maturity':
            if weight > Results.loc['Top Weight, Maturity (g/m^2)',Plant_ID]:
                Results.loc['Top Weight, Maturity (g/m^2)',Plant_ID] = weight
            if Results.loc['Grain Pool, Maturity (g/m^2)',Plant_ID] == 0 and row['Stage Sum Degree Days'] > 55:
                Results.loc['Grain Pool, Maturity (g/m^2)',Plant_ID] = assimilatePool

        #Graphing
        # tillerData.append(dry_matter['N_n'].dropna().values[-1])
        # if row['Stage'] == 'Double Ridge':
        #     if NnPeak.empty:
        #         NnPeak = dry_matter['N_n']
        #     tillerSurvival.append((row['Stage Sum Degree Days'],dry_matter['Proportion Surviving'].dropna()))
        if LAI_z.empty is False:
            peakLAI = max(peakLAI, LAI_z['LAI'].values[-1])
            LAIGraph.append((row['Total Degree Days'],peakLAI))
        
        print("=====================")
        print(row[['Date','Stage','Stage Sum Degree Days']])
        print(dry_matter[['Cohort','#Tillers','Leaf Number','Leaf Active Area']])
        print(LAI_z)

    plant_data['Sum Unaffected Daily Thermal Time'] = plant_data['Daily Unaffected Thermal Time'].cumsum()
    plant_data[['Date','Stage','Total Degree Days','Sum Unaffected Daily Thermal Time']].to_csv(os.path.join(project_path,'Data','Processed','Thermal Time',f'{Plant_ID}.csv'),index=False)
    #raise #just for testing without having to do all 6 plantings
    #plot_tillers(tillerSurvival,NnPeak)
    plot_peakLAI(LAIGraph)

print(Results)
