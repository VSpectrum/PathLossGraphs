from sympy import symbols,Eq,solve,log
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot

from math import pi as pi
import numpy as np
import matplotlib.pyplot as plt

PL, f, hb, hr, d, C = symbols('PL f hb hr d C')
FSPL = 20*log(d, 10)+20*log(f,10)+32.44 #f-GHz, d-m

# f in GHz, d in m, h in m
#M2135_LOS = (20*log((40*pi*d*f/3),10) + 1.5744*log(d,10) - 2.30915 + 0.002*d) * ((2*pi*hb*hr*f*10**9)/(3*10**8)) #+ 40*log(d/((2*pi*hb*hr*f*10**9)/(3*10**8)),10)

#Winner2_NLOS = A*log(d,10) + B + C*log((f/5.0),10) + X
#d in m
Winner2_LOS = 40.0*log(d,10) + 13.47 - 14.0*log(hb,10) - 14.0*log(hr,10) + 6*log((f/5.0),10)

#COST231
#f is in MHz, d in km
COST231_LOS = 42.6 + 26*log(d,10) + 20*log(f,10)

#Ericsson
#f in GHz, h~ in m, d in km
Ericsson = 36.2 + 30.2*log(d,10) + 12*log(hb,10) + 0.1*log(hb,10)*log(d,10) - 3.2*( (log(11.75*hr, 10))**2) + (44.49*log(f,10) - (4.78*(log(f,10))**2) )


#ECC33
#d is km, f is GHz, hb is m, hr is m
ECC33_large = 92.4 + 20*log(d,10) + 20*log(f,10) + 20.41 + 9.83*log(d,10) + 7.894*log(f,10) + 9.56*((log(f,10))**2) - ( (log((hb/200.0),10))*(13.958 + 5.8*((log(d,10))**2)) ) - (0.759*hr - 1.892)

#--------------------------------------------------------------------------------------
FSPL_eval = FSPL.subs([(f,3.5)]).evalf()
#M2135LOS_eval = M2135_LOS.subs([(f,3.5), (hb,35), (hr,1.5)]).evalf()
Winner2LOS_eval = Winner2_LOS.subs([(f,3.5), (hb,35), (hr,1.5)]).evalf()
COST231_LOS_eval = COST231_LOS.subs([(f,3500)]).evalf()
Ericsson_eval = Ericsson.subs([(f,3.5), (hb,35), (hr,1.5)]).evalf()
ECC33_large_eval = ECC33_large.subs([(f,3.5), (hb,35), (hr,1.5)]).evalf()
#--------------------------------------------------------------------------------------
lb_FSPL = lambdify(d, FSPL_eval, modules=['numpy'])
lb_Winner2LOS = lambdify(d, Winner2LOS_eval, modules=['numpy'])
#lb_M2135LOS = lambdify(d, M2135LOS_eval, modules=['numpy'])
lb_COST231_LOS = lambdify(d, COST231_LOS_eval, modules=['numpy'])
lb_Ericsson = lambdify(d, Ericsson_eval, modules=['numpy'])
lb_ECC33_large = lambdify(d, ECC33_large_eval, modules=['numpy'])
#--------------------------------------------------------------------------------------
#d is 0.3 to 5km
t = np.linspace(100, 5000, 10000)
t_km = [x/1000 for x in t]
t_km = np.array(t_km)

t1 = np.linspace(0.1, 5, 10000)
#--------------------------------------------------------------------------------------
plt.plot(t_km, lb_FSPL(t), color = 'black', label='FSPL')
#plt.plot(t_km, lb_M2135LOS(t), color ='deeppink', label='M.2135 LOS')
plt.plot(t_km, lb_Winner2LOS(t), color = 'blue', label='Winner2 LOS')

plt.plot(t1, lb_Ericsson(t1), color = 'lawngreen', label='Ericsson')
plt.plot(t1, lb_ECC33_large(t1), color ='green', label='ECC-33')
plt.plot(t1, lb_COST231_LOS(t1), color ='orange', label='COST 231')

#--------------------------------------------------------------------------------------
plt.title('Path Loss between LOS Urban models')
plt.xlabel('Distance (km)')
plt.ylabel('Path Loss (dB)')

plt.legend(loc='best')
plt.show()