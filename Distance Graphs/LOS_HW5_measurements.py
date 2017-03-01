from sympy import symbols,Eq,solve,log
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot

#from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from scipy import log as slog

import numpy as np
import matplotlib.pyplot as plt

PL, f, hb, hr, d, C = symbols('PL f hb hr d C')

#M2135 Indoor Hotspot www.itu.int/dms_pub/itu-r/opb/rep/R-REP-M.2135-1-2009-PDF-E.pdf pg30
M2135_LOS = 16.9*log(d,10) + 20*log(f,10) + 32.8 #f in GHz, d in m

Lf = symbols('Lf') #0dB for 0 floors penetrated, else 18dB for 1 floor, 26 for 2 floor.
P1238 = 20*log(f,10) + 27*log(d,10) + Lf -28

Winner2_LOS = 18.7*log(d,10) + 46.8 + 20*log((f/5.0),10)

FSPL = 20*log(d, 10)+20*log(f,10)+32.44
#--------------------------------------------------------------------------------------
P1238_eval = P1238.subs([(f,3500), (Lf, 0)]).evalf()
M2135LOS_eval = M2135_LOS.subs([(f,3.5)]).evalf()
Winner2_LOS_eval = Winner2_LOS.subs([(f,3.5)]).evalf()
FSPL_eval = FSPL.subs([(f,3.5)]).evalf()
#--------------------------------------------------------------------------------------
lb_P1238 = lambdify(d, P1238_eval, modules=['numpy'])
lb_M2135LOS = lambdify(d, M2135LOS_eval, modules=['numpy'])
lb_Winner2LOS = lambdify(d, Winner2_LOS_eval, modules=['numpy'])
lb_FSPL = lambdify(d, FSPL_eval, modules=['numpy'])
#--------------------------------------------------------------------------------------
t = [0.4572, 0.762, 1.31064, 1.88976, 0.6096, 1.2192, 1.8288, 2.4384, 3.048, 0.85344, 1.3716, 1.92024, 2.49936, 3.10896, 1.3716, 1.73736, 2.19456, 2.71272, 3.29184, 1.92024, 2.19456, 2.5908, 3.048, 3.56616, 3.048, 3.44424, 3.56616, 4.29768, 3.6576, 3.77952, 3.9624, 4.23672, 4.572, 4.96824, 5.39496, 4.2672, 4.35864, 4.54152, 4.75488, 5.05968, 5.42544, 5.82168, 4.8768, 4.96824, 5.12064, 5.60832, 5.91312, 6.27888, 5.4864, 5.54736, 5.69976, 6.12648, 6.43128, 6.76656, 6.096, 6.15696, 6.27888, 6.67512, 6.94944, 7.28472, 6.7056, 6.76656, 6.88848, 7.25424, 7.4676, 7.80288, 7.3152, 7.37616, 7.4676, 7.62, 7.80288, 8.04672, 8.32104, 7.9248, 7.98576, 8.0772, 8.19912, 8.382, 8.59536, 8.8392, 9.144, 9.47928, 9.81456, 8.62584, 8.96112, 9.32688]
measured_power = [-70.8, -71.5, -71.1, -79.7, -77.9, -78.3, -78.8, -80.6, -81.5, -80.5, -81.3, -81.6, -81.7, -82.3, -71.8, -72.4, -72.9, -74.8, -76.7, -80.4, -81.1, -82.6, -82.8, -83.2, -82.4, -81.4, -80.9, -81.2, -81.8, -82.3, -82.7, -82.6, -83.4, -85.5, -86.3, -86.4, -86.2, -86.4, -86.9, -86.6, -85.7, -86.2, -86.6, -86.7, -86.8, -87.3, -87.4, -87.6, -82.85, -85.25, -86.15, -86.95, -87.65, -87.95, -88.35, -87.35, -86.35, -86.75, -86.95, -87.65, -87.55, -87.65, -87.85, -88.65, -85.5, -85.5, -86.25, -87.25, -87.35, -86.45, -86.35, -86.85, -87.65, -89.15, -89.65, -89.45, -89.65, -89.95, -90.15, -90.35, -93.15, -92.95, -93.25, -93.95, -93.85, -93.55]
#Emitted power is -30dB
D = np.linspace(0.1, 10, 1000)
path_loss_in_NL = [-30-pl for pl in measured_power] #NL - Nortel

#--------------------------------------------------------------------------------------
#PLOTTING CURVE OF BEST-FIT for LOG GRAPH
dictionary = dict(zip(t, path_loss_in_NL))
items = dictionary.items()

x = [a for a,b in items]
y = [b for a,b in items]

x = np.array(x)
y = np.array(y)

def func(x, a, b, c):
    return a*slog(b*x)+c #<--Model for a log curve a*log(b*x)+c -- Best-fitting to this log equation!

popt, pcov = curve_fit(func, x, y)

xs = D
ys = func(xs, *popt)

'''
order = np.argsort(x)
s = UnivariateSpline(x[order], y[order], s=10000, k=5)
xs = np.linspace(min(x), max(x), 1000)
ys = s(xs)
plt.plot(xs, ys,color='black',label="Best Curve")
'''

#--------------------------------------------------------------------------------------
plt.plot(D, lb_FSPL(D), color = 'black', label='FSPL')
plt.plot(D, lb_P1238(D), 'deeppink', label='P.1238')
plt.plot(D, lb_M2135LOS(D), 'lime', linewidth=3, label='M.2135 LOS')
plt.plot(D, lb_Winner2LOS(D), 'blue', label='Winner II LOS')
plt.plot(t, path_loss_in_NL, '*', label='Path Loss in Nortel Lab (Tx=-30dB)')
plt.plot(xs, ys, '--', linewidth=3, color='indigo', label="Best-Fit Curve")
#--------------------------------------------------------------------------------------
plt.title('Measured Indoor Path Loss and Indoor LOS Models for different Distances')
plt.xlabel('Distance (m)')
plt.ylabel('Path Loss (dB)')

plt.legend(loc='best')
plt.show()