from sympy import symbols,Eq,solve,log
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot

from math import pi as pi
import numpy as np
import matplotlib.pyplot as plt

PL, f, hb, hr, d, C = symbols('PL f hb hr d C')

# f in GHz, d in m, h in m
#M2135_LOS = (20*log((40*pi*d*f/3),10) + 1.5744*log(d,10) - 2.30915 + 0.002*d) * ((2*pi*hb*hr*f*10**9)/(3*10**8)) #+ 40*log(d/((2*pi*hb*hr*f*10**9)/(3*10**8)),10)
d0 = 100
c_light = 3*10**8
a, b, c = symbols('a b c')
s_fade = 8.5 #can be between 8.2 to 10.6 dB accounting for shadow fading by trees
SUI_TerrainAB = 20*log((4*pi*d0*(f*10**6)/c_light),10) + 10*(a - b*hb + (c/hb))*log(d/d0,10) + 6*log((f/2000),10) - (10.8*log((hr/2000),10)) + s_fade
SUI_TerrainA_eval = SUI_TerrainAB.subs([(f,3500), (d,1000), (hr,1.5), (a,4.6), (b,0.0075), (c,12.6)]).evalf()
lb_SUI_TerrainA = lambdify(hb, SUI_TerrainA_eval, modules=['numpy'])

#Winner2_NLOS = A*log(d,10) + B + C*log((f/5.0),10) + X
#d in m
Winner2_LOS = 40.0*log(d,10) + 11.65 - 16.2*log(hb,10) - 16.2*log(hr,10) + 3.8*log((f/5.0),10)

#Ericsson
#f in GHz, h~ in m, d in km
Ericsson = 43.2 + 68.93*log(d,10) + 12*log(hb,10) + 0.1*log(hb,10)*log(d,10) - 3.2*( (log(11.75*hr, 10))**2) + (44.49*log(f,10) - (4.78*(log(f,10))**2) )

#ECC33
#d is km, f is GHz, hb is m, hr is m
ECC33_medium = 92.4 + 20*log(d,10) + 20*log(f,10) + 20.41 + 9.83*log(d,10) + 7.894*log(f,10) + 9.56*((log(f,10))**2) - ( (log((hb/200.0),10))*(13.958 + 5.8*((log(d,10))**2)) ) - ( (42.57 + 13.7*log(f,10))*((log(hr,10))-0.585) )

#--------------------------------------------------------------------------------------
Winner2LOS_eval = Winner2_LOS.subs([(f,3.5), (d,1000), (hr,1.5)]).evalf()
Ericsson_eval = Ericsson.subs([(f,3.5), (d,1), (hr,1.5)]).evalf()
ECC33_medium_eval = ECC33_medium.subs([(f,3.5), (d,1), (hr,1.5)]).evalf()
#--------------------------------------------------------------------------------------
lb_Winner2LOS = lambdify(hb, Winner2LOS_eval, modules=['numpy'])
lb_Ericsson = lambdify(hb, Ericsson_eval, modules=['numpy'])
lb_ECC33_medium = lambdify(hb, ECC33_medium_eval, modules=['numpy'])
#--------------------------------------------------------------------------------------
t = np.linspace(10, 80, 1000)
#--------------------------------------------------------------------------------------
plt.plot(t, lb_Winner2LOS(t), color = 'blue', label='Winner2 LOS')
plt.plot(t, lb_Ericsson(t), color = 'red', label='Ericsson')
plt.plot(t, lb_ECC33_medium(t), color ='green', label='ECC-33')
plt.plot(t, lb_SUI_TerrainA(t), '--', color = 'darkorchid', label='SUI Terrain A')
#--------------------------------------------------------------------------------------
plt.title('Path Loss between LOS Suburban models')
plt.xlabel('Base Station Antenna Height (m)')
plt.ylabel('Path Loss (dB)')

plt.legend(loc='best')
plt.show()