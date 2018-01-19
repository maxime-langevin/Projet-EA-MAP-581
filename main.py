
import numpy as np
import matplotlib.pyplot as plt
from scipy . integrate import odeint

nbre_points = 3001
t_tab = np.linspace(0 , 500 , nbre_points)

A = 1
P = 0.1
alpha_E = 0.05 
alpha_I = 0.05
w = 0.05
phi = 0

E0 = 0.2
I0 = 0.2

wEE = 27
wEI = 50
wIE = 15
wII = 0

def SE(XE, XI, P, E0, wEE, wEI, wII, wIE):
    exp_value = np.exp(-wEE*XE + wIE*XI - P)
    D = 1 + ( 1 / E0 - 2 ) * exp_value
    return 1 / ( 1 + D )

def SI(XE, XI, I0, wEE, wEI, wII, wIE):
    exp_value = np.exp(-wEI*XE + wII*XI)
    D = 1 + ( 1 / I0 - 2 ) * exp_value
    return 1 / ( 1 + D )

def XE_prime(XE, XI, P, E0, wEE, wEI, wII, wIE):
    return -E0 - XE + (1 - E0 -XE) * SE(XE, XI, P, E0, wEE, wEI, wII, wIE)

def XI_prime(XE, XI, A, I0, wEE, wEI, wII, wIE):
    return 1/A * (-I0 - XI + (1 - I0 -XI) * SI(XE, XI, I0, wEE, wEI, wII, wIE))

def X_prime(X, t, P, A, E0, I0, wEE, wEI, wII, wIE, alpha_E, alpha_I, w, phi):
    XE, XI = X
    E0_t = E0 + alpha_E * np.cos(w * t)
    I0_t = I0 + alpha_I * np.cos(w * t + phi)
    rep = [XE_prime(XE, XI, P, E0_t, wEE, wEI, wII, wIE), XI_prime(XE, XI, A, I0_t, wEE, wEI, wII, wIE)]
    return rep

print(t_tab)
sol = odeint(X_prime, [0, 0], t_tab, args = (P, A, E0, I0, wEE, wEI, wII, wIE, alpha_E, alpha_I, w, phi))

plt.figure()
plt.plot()
plt.plot(t_tab, sol[:, 0])
plt.plot(t_tab, E0 + alpha_E * np.cos(w * t_tab))
#plt.plot(t_tab, sol[:, 1])
plt.savefig("Fig1")
plt.show()
