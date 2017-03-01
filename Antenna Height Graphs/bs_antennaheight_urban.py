from sympy import symbols,Eq,solve,log
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot

from math import pi as pi
import numpy as np
import matplotlib.pyplot as plt

PL, f, hb, hr, d, C = symbols('PL f hb hr d C')

#Winner2_NLOS = A*log(d,10) + B + C*log((f/5.0),10) + X
#d in m
Winner2_LOS = 40.0*log(d,10) + 13.47 - 14.0*log(hb,10) - 14.0*log(hr,10) + 6*log((f/5.0),10)

#Ericsson
#f in GHz, h~ in m, d in km
Ericsson = 36.2 + 30.2*log(d,10) + 12*log(hb,10) + 0.1*log(hb,10)*log(d,10) - 3.2*( (log(11.75*hr, 10))**2) + (44.49*log(f,10) - (4.78*(log(f,10))**2) )

#ECC33
#d is km, f is GHz, hb is m, hr is m
ECC33_large = 92.4 + 20*log(d,10) + 20*log(f,10) + 20.41 + 9.83*log(d,10) + 7.894*log(f,10) + 9.56*((log(f,10))**2) - ( (log((hb/200.0),10))*(13.958 + 5.8*((log(d,10))**2)) ) - (0.759*hr - 1.892)

#--------------------------------------------------------------------------------------
Winner2LOS_eval = Winner2_LOS.subs([(f,3.5), (d,1000), (hr,1.5)]).evalf()
Ericsson_eval = Ericsson.subs([(f,3.5), (d,1), (hr,1.5)]).evalf()
ECC33_large_eval = ECC33_large.subs([(f,3.5), (d,1), (hr,1.5)]).evalf()
#--------------------------------------------------------------------------------------
lb_Winner2LOS = lambdify(hb, Winner2LOS_eval, modules=['numpy'])
lb_Ericsson = lambdify(hb, Ericsson_eval, modules=['numpy'])
lb_ECC33_large = lambdify(hb, ECC33_large_eval, modules=['numpy'])
#--------------------------------------------------------------------------------------
#d is 0.3 to 5km
t = np.linspace(10, 80, 1000)
#--------------------------------------------------------------------------------------
plt.plot(t, lb_Winner2LOS(t), color = 'blue', label='Winner2 LOS')
plt.plot(t, lb_Ericsson(t), color = 'red', label='Ericsson')
plt.plot(t, lb_ECC33_large(t), color ='green', label='ECC-33')
#--------------------------------------------------------------------------------------
plt.title('Path Loss between LOS Urban models')
plt.xlabel('Base Station Antenna Height (m)')
plt.ylabel('Path Loss (dB)')

plt.legend(loc='best')
plt.show()