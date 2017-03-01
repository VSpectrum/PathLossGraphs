from sympy import symbols,Eq,solve,log
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot

import numpy as np
import matplotlib.pyplot as plt

PL, f, hb, hr, d, C = symbols('PL f hb hr d C')
FSPL = 20*log(d, 10)+20*log(f,10)+32.44 #f-GHz, d-m
Hata = (46.3+33.9*log(f,10)-13.82*log(hb,10)-(hr*(1.1*log(f,10)-0.7)-(1.56*log(f,10)-0.8))+log(d,10)*(44.9-6.55*log(hb,10))+C)

#ECC33
#d is km, f is Ghz, hb is m, hr is m
ECC33_large = 92.4 + 20*log(d,10) + 20*log(f,10) + 20.41 + 9.83*log(d,10) + 7.894*log(f,10) + 9.56*((log(f,10))**2) - ( (log((hb/200.0),10))*(13.958 + 5.8*((log(d,10))**2)) ) - (0.759*hr - 1.892)
ECC33_medium = 92.4 + 20*log(d,10) + 20*log(f,10) + 20.41 + 9.83*log(d,10) + 7.894*log(f,10) + 9.56*((log(f,10))**2) - ( (log((hb/200.0),10))*(13.958 + 5.8*((log(d,10))**2)) ) - ( (42.57 + 13.7*log(f,10))*((log(hr,10))-0.585) )

#Ericsson
a0, a1, a2, a3 = symbols('a0 a1 a2 a3') #a0 = 36.2, a1=30.2, a2 = 12, a3 = 0.1 <- urban
Ericsson = a0 + a1*log(d,10) + a2*log(hb,10) + a3*log(hb,10)*log(d,10) - 3.2*( (log(11.75*hr, 10))**2) + (44.49*log(f,10) - (4.78*(log(f,10))**2) )

#M2135 Indoor Hotspot www.itu.int/dms_pub/itu-r/opb/rep/R-REP-M.2135-1-2009-PDF-E.pdf pg30
M2135_LOS = 16.9*log(d,10) + 20*log(f, 10) + 32.8 #f in GHz, d in m
M2135_NLOS = 43.3*log(d,10) + 20*log(f,10) + 11.5

#WINNER 2? www.raymaps.com/index.php/winner-ii-path-loss-model/
A, B, C, X = symbols('A B C X') #f in GHz, d in m
Winner2_NLOS = A*log(d,10) + B + C*log((f/5.0),10) + X
Winner2_LOS = 18.7*log(d,10) + 46.8 + 20*log((f/5.0),10)

#ITU Indoor/P1238- www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.1238-7-201202-S!!PDF-E.pdf
#N is 27 for the 3.5GHz range in N*log(d,10) , f in MHz, d in m, LF in dB
Lf = symbols('Lf') #0dB for 0 floors penetrated, else 18dB for 1 floor, 26 for 2 floor.
P1238 = 20*log(f,10) + 27*log(d,10) + Lf -28
#--------------------------------------------------------------------------------------
FSPL_eval = FSPL.subs([(f,3.5)]).evalf()
P1238_eval = P1238.subs([(f,3500), (Lf, 0)]).evalf()
M2135LOS_eval = M2135_LOS.subs([(f,3.5)]).evalf()
M2135NLOS_eval = M2135_NLOS.subs([(f,3.5)]).evalf()
Winner2LOS_eval = Winner2_LOS.subs([(f,3.5)]).evalf()
#--------------------------------------------------------------------------------------
lb_FSPL = lambdify(d, FSPL_eval, modules=['numpy'])
lb_P1238 = lambdify(d, P1238_eval, modules=['numpy'])
lb_Winner2LOS = lambdify(d, Winner2LOS_eval, modules=['numpy'])
lb_M2135NLOS = lambdify(d, M2135NLOS_eval, modules=['numpy'])
lb_M2135LOS = lambdify(d, M2135LOS_eval, modules=['numpy'])
#--------------------------------------------------------------------------------------
t = np.linspace(3, 100, 10000)
#--------------------------------------------------------------------------------------
plt.plot(t, lb_FSPL(t), 'black', label='FSPL')
plt.plot(t, lb_P1238(t), 'deeppink', label='P1238')
plt.plot(t, lb_Winner2LOS(t), 'blue', label='Winner2 LOS')
plt.plot(t, lb_M2135NLOS(t), 'lawngreen', label='M2135 NLOS')
plt.plot(t, lb_M2135LOS(t), 'lime', label='M2135 LOS')
#--------------------------------------------------------------------------------------
plt.title('Path Loss vs Distance with Different Propagation Models')
plt.xlabel('Distance (m)')
plt.ylabel('Path Loss (dB)')

plt.legend(loc='best')
plt.show()