from math import pi, sin, cos, asin
import matplotlib.pyplot as plt
import os

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), r'../..'))

#Calculate Declination Angle
def calc_dec_alt(julian_day):
    tilt_of_earth = 23.44 * (pi / 180) #Converts to radians
    radians_per_day = 2*pi * (1 / 365.25)
    declination = -tilt_of_earth * cos(radians_per_day * (julian_day + 10))
    
    return declination

def calc_dec_paper(julian_day):
    theta_1 = 2*pi*(julian_day - 80) / 365
    theta_2 = 0.0335*(sin(2*pi*julian_day) - sin(2*pi*80)) / 365
    theta = theta_1 + theta_2

    declination = asin(0.3978*sin(theta))
    
    return declination

julian_day = [jday for jday in range(366)]

# Compute difference between equations
difference = [calc_dec_alt(jday) - calc_dec_paper(jday) for jday in julian_day]

# Find max and min differences
max_diff = max(difference)
min_diff = min(difference)
max_day = julian_day[difference.index(max_diff)]
min_day = julian_day[difference.index(min_diff)]

# Plot lines
plt.plot(julian_day, [calc_dec_alt(jday) for jday in julian_day], label='Alternate Equation')
plt.plot(julian_day, [calc_dec_paper(jday) for jday in julian_day], label='Original Equation')
plt.plot(julian_day, difference, label='Alternate Equation - Original Equation')

# Annotate max and min points
plt.annotate(f'Max: {max_diff:.4f}', xy=(max_day, max_diff), xytext=(max_day-35, -0.1),
             arrowprops=dict(facecolor='green', arrowstyle='->'), fontsize=11)
plt.annotate(f'Min: {min_diff:.4f}', xy=(min_day, min_diff), xytext=(min_day-35, 0.1),
             arrowprops=dict(facecolor='red', arrowstyle='->'), fontsize=11)

# Plot formatting
plt.xlabel('Day of Year')
plt.ylabel('Declination angle (radians)')
plt.title('Declination Angle Equation Comparison')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(project_path, 'Data/Graphs', 'DeclinationEquationComparison.png'))
plt.show()
plt.close()