from sympy import symbols,Eq,solve,log
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot
from math import pi
from scipy.interpolate import InterpolatedUnivariateSpline as UnivariateSpline

import numpy as np
import matplotlib.pyplot as plt

#SUI Terrain B, COST-231, Winner 2, M.2135

PL, f, hb, hr, d, C = symbols('PL f hb hr d C')

#M2135 Indoor Hotspot www.itu.int/dms_pub/itu-r/opb/rep/R-REP-M.2135-1-2009-PDF-E.pdf pg30
M2135_LOS = 16.9*log(d,10) + 20*log(f,10) + 32.8 #f in GHz, d in m
Winner2_LOS = 18.7*log(d,10) + 46.8 + 20*log((f/5.0),10)
FSPL = 20*log(d, 10)+20*log(f,10)+32.44

d0=100
c_light=3*10**8
a, b, c = symbols('a b c')
s_fade = 8.5 
SUI_TerrainAB = 20*log((4*pi*d0*(f*10**6)/c_light),10) + 10*(a - b*hb + (c/hb))*log(d/d0,10) + 6*log((f/2000),10) + (-10.8*log((hr/2000),10)) + s_fade
COST231_LOS = 42.6 + 26*log(d,10) + 20*log(f,10)
#--------------------------------------------------------------------------------------

M2135LOS_eval = M2135_LOS.subs([(f,3.5)]).evalf()
Winner2_LOS_eval = Winner2_LOS.subs([(f,3.5)]).evalf()
FSPL_eval = FSPL.subs([(f,3.5)]).evalf()

SUI_TerrainB_eval = SUI_TerrainAB.subs([(f,3500), (hb,35), (hr,1.5), (a,4), (b,0.0065), (c,17.1)]).evalf()
COST231_LOS_eval = COST231_LOS.subs([(f,3500)]).evalf()
#--------------------------------------------------------------------------------------
lb_COST231_LOS = lambdify(d, COST231_LOS_eval, modules=['numpy'])
lb_SUI_TerrainB = lambdify(d, SUI_TerrainB_eval, modules=['numpy'])
lb_M2135LOS = lambdify(d, M2135LOS_eval, modules=['numpy'])
lb_Winner2LOS = lambdify(d, Winner2_LOS_eval, modules=['numpy'])
lb_FSPL = lambdify(d, FSPL_eval, modules=['numpy'])
#--------------------------------------------------------------------------------------

t = np.linspace(50, 100, 100)
t_km = [x/1000 for x in t]
t_km = np.array(t_km)
#--------------------------------------------------------------------------------------
t1 = np.linspace(0.1, 0.2, 100)
#--------------------------------------------------------------------------------------
plt.plot(t, lb_COST231_LOS(t1), color ='orange', label='COST 231')
plt.plot(t, lb_SUI_TerrainB(t), color = 'red', label='SUI Terrain B')
plt.plot(t, lb_FSPL(t), color = 'black', label='FSPL')
plt.plot(t, lb_M2135LOS(t), 'lime', label='M2135 LOS')
plt.plot(t, lb_Winner2LOS(t), 'blue', label='Winner2 LOS')
#--------------------------------------------------------------------------------------
plt.title('Indoor Models and Outdoor Models')
plt.xlabel('Distance (m)')
plt.ylabel('Path Loss (dB)')

plt.legend(loc='best')
plt.show()