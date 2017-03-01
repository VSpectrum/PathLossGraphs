from sympy import symbols,Eq,solve,log
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot

from math import pi as pi
import numpy as np
import matplotlib.pyplot as plt

PL, f, hb, hr, d, C = symbols('PL f hb hr d C')
FSPL = 20*log(d, 10)+20*log(f,10)+32.44 #f-GHz, d-m

# f in GHz, d in m, h in m
M2135_LOS = (20*log((40*pi*d*f/3),10) + 0.477*log(d,10) - 0.7 + 0.002*d*0.699) * ((2*pi*hb*hr*f)/(3)) + 40*log(d/((2*pi*hb*hr*f)/(3)),10)

#SUI f in MHz, d in m
'''
Parameter 	Ter A 	Ter B 	Ter C 
a 			4.6 	4.0 	3.6 
b(m-1) 		0.0075 	0.0065 	0.005 
c(m) 		12.6 	17.1 	20 
'''
c_light = 3*10**8
d0 = 100
a, b, c = symbols('a b c')
s_fade = 8.5 #can be between 8.2 to 10.6 dB accounting for shadow fading by trees
		#																						Xf				Xh
SUI_TerrainC = 20*log((4*pi*d0*(f*10**6)/c_light),10) + 10*(a - b*hb + (c/hb))*log(d/d0,10) + 6*log((f/2000),10) + (-20*log((hr/2000),10)) + s_fade
SUI_TerrainAB = 20*log((4*pi*d0*(f*10**6)/c_light),10) + 10*(a - b*hb + (c/hb))*log(d/d0,10) + 6*log((f/2000),10) + (-10.8*log((hr/2000),10)) + s_fade

#Winner2_NLOS = A*log(d,10) + B + C*log((f/5.0),10) + X  #d in m
Winner2_LOS = 40.0*log(d,10) + 10.5 - 18.5*log(hb,10) - 18.5*log(hr,10) + 1.5*log((f/5.0),10)

#Ericsson  #f in GHz, h~ in m, d in km
Ericsson = 45.95 + 100.6*log(d,10) + 12*log(hb,10) + 0.1*log(hb,10)*log(d,10) - 3.2*( (log(11.75*hr, 10))**2) + (44.49*log(f,10) - (4.78*(log(f,10))**2) )

#--------------------------------------------------------------------------------------
FSPL_eval = FSPL.subs([(f,3.5)]).evalf()
M2135LOS_eval = M2135_LOS.subs([(f,3.5), (hb,35), (hr,1.5)]).evalf()
Winner2LOS_eval = Winner2_LOS.subs([(f,3.5), (hb,35), (hr,1.5)]).evalf()
Ericsson_eval = Ericsson.subs([(f,3.5), (hb,35), (hr,1.5)]).evalf()
SUI_TerrainA_eval = SUI_TerrainAB.subs([(f,3500), (hb,35), (hr,1.5), (a,4.6), (b,0.0075), (c,12.6)]).evalf()
SUI_TerrainB_eval = SUI_TerrainAB.subs([(f,3500), (hb,35), (hr,1.5), (a,4), (b,0.0065), (c,17.1)]).evalf()
SUI_TerrainC_eval = SUI_TerrainC.subs([(f,3500), (hb,35), (hr,1.5), (a,3.6), (b,0.005), (c,20)]).evalf()
#--------------------------------------------------------------------------------------
lb_FSPL = lambdify(d, FSPL_eval, modules=['numpy'])
lb_Winner2LOS = lambdify(d, Winner2LOS_eval, modules=['numpy'])
lb_M2135LOS = lambdify(d, M2135LOS_eval, modules=['numpy'])
lb_Ericsson = lambdify(d, Ericsson_eval, modules=['numpy'])
lb_SUI_TerrainA = lambdify(d, SUI_TerrainA_eval, modules=['numpy'])
lb_SUI_TerrainB = lambdify(d, SUI_TerrainB_eval, modules=['numpy'])
lb_SUI_TerrainC = lambdify(d, SUI_TerrainC_eval, modules=['numpy'])

#--------------------------------------------------------------------------------------
#d is 0.3 to 5km
t = np.linspace(100, 5000, 90)
t_km = [x/1000 for x in t]
t_km = np.array(t_km)

t1 = np.linspace(0.1, 5, 90)
#--------------------------------------------------------------------------------------
plt.plot(t_km, lb_FSPL(t), color = 'black', label='FSPL')
#plt.plot(t_km, lb_M2135LOS(t), color ='deeppink', label='M.2135 LOS') #Oddity with equation provided in the M.2135 ITU-R manual producing incorrect graph
plt.plot(t_km, lb_Winner2LOS(t), color = 'blue', label='Winner2 LOS')
plt.plot(t1, lb_Ericsson(t1), color = 'red', label='Ericsson')

plt.plot(t_km, lb_SUI_TerrainA(t), '--', color = 'blue', label='SUI Terrain A')
plt.plot(t_km, lb_SUI_TerrainB(t), '--', color = 'red', label='SUI Terrain B')
plt.plot(t_km, lb_SUI_TerrainC(t), '--', color = 'green', label='SUI Terrain C')
#--------------------------------------------------------------------------------------
#plt.grid(b=True, which='major', color='gray', linestyle='--')
plt.title('Path Loss between LOS Rural models')
plt.xlabel('Distance (km)')
plt.ylabel('Path Loss (dB)')

plt.legend(loc='best')
plt.show()