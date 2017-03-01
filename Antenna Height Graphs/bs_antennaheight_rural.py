from sympy import symbols,Eq,solve,log
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot

from math import pi as pi
import numpy as np
import matplotlib.pyplot as plt

PL, f, hb, hr, d, C = symbols('PL f hb hr d C')

#SUI f in MHz, d in m
c_light = 3*10**8
d0 = 100
a, b, c = symbols('a b c')
s_fade = 8.5 #can be between 8.2 to 10.6 dB accounting for shadow fading by trees
SUI_TerrainC = 20*log((4*pi*d0*(f*10**6)/c_light),10) + 10*(a - b*hb + (c/hb))*log(d/d0,10) + 6*log((f/2000),10) + (-20*log((hr/2000),10)) + s_fade
SUI_TerrainAB = 20*log((4*pi*d0*(f*10**6)/c_light),10) + 10*(a - b*hb + (c/hb))*log(d/d0,10) + 6*log((f/2000),10) + (-10.8*log((hr/2000),10)) + s_fade

#Winner2_NLOS = A*log(d,10) + B + C*log((f/5.0),10) + X  #d in m
Winner2_LOS = 40.0*log(d,10) + 10.5 - 18.5*log(hb,10) - 18.5*log(hr,10) + 1.5*log((f/5.0),10)

#Ericsson  #f in GHz, h~ in m, d in km
Ericsson = 45.95 + 100.6*log(d,10) + 12*log(hb,10) + 0.1*log(hb,10)*log(d,10) - 3.2*( (log(11.75*hr, 10))**2) + (44.49*log(f,10) - (4.78*(log(f,10))**2) )

#--------------------------------------------------------------------------------------
Winner2LOS_eval = Winner2_LOS.subs([(f,3.5), (d,1000), (hr,1.5)]).evalf()
Ericsson_eval = Ericsson.subs([(f,3.5), (d,1), (hr,1.5)]).evalf()
SUI_TerrainA_eval = SUI_TerrainAB.subs([(f,3500), (d,1000), (hr,1.5), (a,4.6), (b,0.0075), (c,12.6)]).evalf()
SUI_TerrainB_eval = SUI_TerrainAB.subs([(f,3500), (d,1000), (hr,1.5), (a,4), (b,0.0065), (c,17.1)]).evalf()
SUI_TerrainC_eval = SUI_TerrainC.subs([(f,3500), (d,1000), (hr,1.5), (a,3.6), (b,0.005), (c,20)]).evalf()
#--------------------------------------------------------------------------------------
lb_Winner2LOS = lambdify(hb, Winner2LOS_eval, modules=['numpy'])
lb_Ericsson = lambdify(hb, Ericsson_eval, modules=['numpy'])
lb_SUI_TerrainA = lambdify(hb, SUI_TerrainA_eval, modules=['numpy'])
lb_SUI_TerrainB = lambdify(hb, SUI_TerrainB_eval, modules=['numpy'])
lb_SUI_TerrainC = lambdify(hb, SUI_TerrainC_eval, modules=['numpy'])
#--------------------------------------------------------------------------------------
t = np.linspace(10, 80, 1000)
#--------------------------------------------------------------------------------------
plt.plot(t, lb_Winner2LOS(t), color = 'blue', label='Winner2 LOS')
plt.plot(t, lb_Ericsson(t), color = 'red', label='Ericsson')

plt.plot(t, lb_SUI_TerrainA(t), '--', color = 'blue', label='SUI Terrain A')
plt.plot(t, lb_SUI_TerrainB(t), '--', color = 'red', label='SUI Terrain B')
plt.plot(t, lb_SUI_TerrainC(t), '--', color = 'green', label='SUI Terrain C')
#--------------------------------------------------------------------------------------
#plt.grid(b=True, which='major', color='gray', linestyle='--')
plt.title('Path Loss between LOS Rural models')
plt.xlabel('Base Station Antenna Height (m)')
plt.ylabel('Path Loss (dB)')

plt.legend(loc='best')
plt.show()