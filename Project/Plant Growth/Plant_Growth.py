from math import acos, asin, cos, exp, pi, sin, sqrt, tan
import datetime

########   Phenological Development Submodel   #########

###Thermal Time
#Manually set variables for testing
T_min = max(int(input("T_min (Celcius): ")),0)
T_max = max(int(input("T_max (Celcius): ")),0)
T_base = max(int(input("T_base (Celcius): ")),0)

#Constants given in paper
T_opt = 26
TD_max = 37

#Formula to calculate Thermal Time (measured in degree days (Celsius))
#Degree days are a measure of how cold or warm a location is
#calculated by comparing the average daily temperature to a standard temperature
T_t = 0
for r in range(1,9):
    #print("iteraion:", r)
    f_r = (1/2)*(1+cos((90/8)*(2*r-1)))
    #print("f_r =", f_r)
    T_H = max(T_min+f_r*(T_max-T_min),0) #Degree Celcius
    #print("T_H =", T_H)
    if T_H < T_opt:
        T_t += T_H-T_base
        #print("T_t increased by:", T_H-T_base)
    elif T_H == T_opt:
        T_t += T_opt-T_base
        #print("T_t increased by:", T_opt-T_base)
    elif T_H < TD_max:
        T_t += (T_opt-T_base)*(TD_max-T_H)/(TD_max-T_opt)
        #print("T_t increased by:", (T_opt-T_base)*(TD_max-T_H)/(TD_max-T_opt))
    else:
        T_t += 0
        #print("T_t increased by: 0")
T_t = (1/8)*T_t #Degree Celcius Days
#print("Thermal Time =", T_t)

###Photoperiod Effects
#Constants given in paper
P_opt = 20
P_Base = 0 #P_base = 0 from emergence to double ridge, but P_Base = 7 from double ridge to anthesis

#Calculating Julian Day (days so far in the current year)
Date_str = input("What is the date? (DD/MM/YYYY)")
Date_dt = datetime.datetime.strptime(Date_str, '%d/%m/%Y')
Date_tt = Date_dt.timetuple()
Jday = Date_tt.tm_yday

#Calculating Declension Angle
theta1 = 2*pi*(Jday-80)/365 #Radians
theta2 = 0.0335*(sin(2*pi*Jday)-sin(2*pi*80)) #Radians
theta = theta1 + theta2 #Radians
Dec = asin(0.3978*sin(theta)) #Radians

Lat = int(input("Latitude: "))
D = -0.10453/(cos(Lat)*cos(Dec))
P_R = acos(D-tan(Lat)*tan(Dec)) #Radians
P_H = 24*P_R/pi #Hours
FP = (P_H-P_Base)/(P_opt-P_Base)


###Vernalization
VDD = 0
V_sat = 33
V_base = 8


n = int(input("number of days from germination: "))
for i in range(1,n+1):
    V_eff = 0
    for r in range(1,9):
        #print("iteraion:", r)
        f_r = (1/2)*(1+cos((90/8)*(2*r-1)))
        #print("f_r =", f_r)
        T_H = max(T_min+f_r*(T_max-T_min),0)
        #print("T_H =", T_H)
        if T_H >= -4 and T_H<3:
            V_eff += (T_H+4)/7
        elif T_H >= 3 and T_H < 10:
            V_eff += 1
        elif T_H >= 10 and T_H < 17:
            V_eff += (17-T_H)/7
    VDD += V_eff #Vernal Days

FV = (VDD-V_base)/(V_sat-V_base)
########################################################

########### Tiller and leaf growth Submodel ############

########################################################

###########   Root Growth Submodel #####################
T = int(input("daily mean temperature"))

TR = 0.2 + 0.12*T

########################################################

#### Light Interception and Photosynthesis Submodel ####

### Radiation Interception
Q_p = int(input("Photosynthetically active radiation (PAR) at the top of the canopy: "))
k = int(input("extinction coefficient: "))
m = int(input("leaf transmission coefficient for PAR: "))
LAI_z = int(input("Green leaf area index at level z (z=0 at canopy): "))

Q_p_z = ((Q_p*k)/(1-m))*exp(-k*LAI_z) #W/m^2

###Phtosynthesis
r_a = int(input("crop boundary layer resistance: "))
r_m = int(input("mesophyll resistance: "))
D = int(input("vapour pressure deficit (kPa): "))
r_s = 1.56*75*(1+100/Q_p_z)*(1-0.3*D)
r_p = r_a + r_s + r_m

C_a = int(input("ambient CO2 concerntration (mg/m^3): "))
alpha = int(input("photosynthetic efficieny (mg/J): "))
P_max = 0.995*(C_a/r_p)

a1 = 0.995
b1 = -(P_max+alpha*Q_p_z)
c1 = alpha*Q_p_z*P_max
P_g1 = (-b1+sqrt(b1**2-4*a1*c1))/(2*a1) #mgCO2/m^2 per sec

###Temperature Correction of maximum photosynthetic rate
T_k = int(input("Leaf Temperature in Kelvin: "))
DeltaH_1 = int(input(" Activation energies for the electron transport system (cal/mol): "))
DeltaH_2 = int(input(" Denaturation energies for the electron transport system (cal/mol): "))
R = 1.987 #gas constant in Cal/K per mol
DeltaS = int(input("entropy change on denaturation of the electron transport system (cal/K per mol): "))

P_m = (0.044*6*(10**9)*T_k*exp(-DeltaH_1/(R*T_k)))/(1+exp(-DeltaH_2/(R*T_k))*exp(DeltaS/R))

a2 = 0.995*(1/P_max)*((1/(alpha*Q_p_z))+(1/P_m))
b2 = -1*((1/(alpha*Q_p_z))+(1/P_m)+1/P_max)
c2 = 1
P_g2 = (-b2+sqrt(b2**2-4*a2*c2))/(2*a2)

###Respiration
alpha2 = int(input("growth rate coefficient: "))
H = int(input("Number of daylight hours: "))
W = int(input("crop weight (g/m^3): "))
b = int(input("maintenance respiration coefficient: "))
P_g_CO2 = 0
for h in range(0,H+1): 
    P_g_CO2 += P_g*(h)+W*b*(2**(0.05*(T_max+T_min)))

R2 = 0.65*alpha2*P_g_CO2 #g/m^2 per day
########################################################

## Dry-matter partitioning and grain growth Submodel ##
P_n = P_g_CH20 - R2 #g/m_2 per day (seems like m_2 is a typo of m^2)

G_max = 0.045*(T_max+T_min)/2 + 0.4 #mg/grain per day
########################################################
