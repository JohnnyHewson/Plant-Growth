from math import cos, pi, sin, asin, acos, tan
import datetime

lat = 77
date = datetime.date(2025,1,1)

Jday = date.timetuple().tm_yday
theta1 = 2 * pi * (Jday - 80) / 365
theta2 = 0.0335 * (sin(2 * pi * Jday) - sin(2 * pi * 80)) / 365
theta = theta1 + theta2
Dec = asin(0.3978 * sin(theta))
D = -0.10453 / (cos(lat * pi / 180) * cos(Dec))
