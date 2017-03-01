from sympy import symbols,Eq,solve,log
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot

import numpy as np
import matplotlib.pyplot as plt

def Winner2_NLOS_X(numwalls):
	return 5.0*(numwalls-1) #12*(numwalls-1) for heavy walls

def Winner2_NLOS_eval(Winner2_NLOS, numwalls):
	return Winner2_NLOS.subs([(f,3.5),(X, Winner2_NLOS_X(numwalls))]).evalf()

PL, f, hb, hr, d, C = symbols('PL f hb hr d C')
FSPL = 20*log(d, 10)+20*log(f,10)+32.44 #f-GHz, d-m

#M2135 Indoor Hotspot www.itu.int/dms_pub/itu-r/opb/rep/R-REP-M.2135-1-2009-PDF-E.pdf pg30
M2135_LOS = 16.9*log(d,10) + 20*log(f, 10) + 32.8 #f in GHz, d in m
M2135_NLOS = 43.3*log(d,10) + 20*log(f,10) + 11.5

#WINNER 2? www.raymaps.com/index.php/winner-ii-path-loss-model/
A, B, C, X = symbols('A B C X') #f in GHz, d in m

#Winner2_NLOS = A*log(d,10) + B + C*log((f/5.0),10) + X
Winner2_NLOS = 36.8*log(d,10) + 43.8 + 20*log((f/5.0),10) + X
Winner2_LOS = 18.7*log(d,10) + 46.8 + 20*log((f/5.0),10)

#ITU Indoor/P1238- www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.1238-7-201202-S!!PDF-E.pdf
#N is 27 for the 3.5GHz range in N*log(d,10) , f in MHz, d in m, LF in dB
Lf = symbols('Lf') #0dB for 0 floors penetrated, else 18dB for 1 floor, 26 for 2 floor.
P1238 = 20*log(f,10) + 27*log(d,10) + Lf -28

FSPL_eval = FSPL.subs([(f,3.5)]).evalf()

M2135LOS_eval = M2135_LOS.subs([(f,3.5)]).evalf()
M2135NLOS_eval = M2135_NLOS.subs([(f,3.5)]).evalf()

#Winner2_NLOS_eval = Winner2_NLOS.subs([(f,3.5),(X, Winner2_NLOS_X(2))]).evalf()
Winner2LOS_eval = Winner2_LOS.subs([(f,3.5)]).evalf()

lb_FSPL = lambdify(d, FSPL_eval, modules=['numpy'])

lb_Winner2NLOS2 = lambdify(d, Winner2_NLOS_eval(Winner2_NLOS, 2), modules=['numpy'])
lb_Winner2NLOS4 = lambdify(d, Winner2_NLOS_eval(Winner2_NLOS, 4), modules=['numpy'])

lb_Winner2LOS = lambdify(d, Winner2LOS_eval, modules=['numpy'])

lb_M2135NLOS = lambdify(d, M2135NLOS_eval, modules=['numpy'])
lb_M2135LOS = lambdify(d, M2135LOS_eval, modules=['numpy'])

t = np.linspace(10, 100, 10000)

plt.plot(t, lb_FSPL(t), color = 'black', label='FSPL')
plt.plot(t, lb_Winner2NLOS2(t), '--', color ='deeppink', label='Winner2 NLOS (2 Walls)')
plt.plot(t, lb_Winner2LOS(t), color = 'blue', label='Winner2 LOS')
plt.plot(t, lb_M2135NLOS(t), '--', color = 'lawngreen', label='M2135 NLOS')
plt.plot(t, lb_M2135LOS(t), color ='lime', label='M2135 LOS')

plt.title('Path Loss between LOS and NLOS Indoor models')
plt.xlabel('Distance (m)')
plt.ylabel('Path Loss (dB)')

plt.legend(loc='best')
plt.show()