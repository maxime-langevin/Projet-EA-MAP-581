
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 17:07:56 2018

@author: maxime
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy . integrate import odeint

nbre_points = 3001
t_tab = np.linspace(0 , 500 , nbre_points)

A = 1

P_sin=0.5
T_sinusoidale=10
phi=np.pi*2
offset=-0.1

Pmax_creneau = 0.5
Pmin_creneau=0
T_creneau=50

P_const=0.5

T=250
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

def X_prime(X, t, A, E0, I0, wEE, wEI, wII, wIE, alpha_E, alpha_I, w, phi,commande):
    XE, XI = X
    E0_t = E0 + alpha_E * np.cos(w * t)
    I0_t = I0 + alpha_I * np.cos(w * t + phi)
    if commande=="creneau":
        P=P_creneau(t,Pmin_creneau,Pmax_creneau,T)
    elif commande=="constante":
        P=P_constante(t,P_const)
    elif commande=="sinusoidale":
        P=P_sinusoidale(t,P_sin,T,phi,offset)
    
    
    rep = [XE_prime(XE, XI, P, E0_t, wEE, wEI, wII, wIE), XI_prime(XE, XI, A, I0_t, wEE, wEI, wII, wIE)]
    return rep

'''Commande Créneau
'''


def P_creneau(t,Pmin,Pmax,T):
    if t%(2*T)>T:
        return Pmin
    else:
        return Pmax
   
        
def P_creneau_tab(t_tab,Pmin,Pmax,T):
    return np.array([Pmax if t%(2*T)>T else Pmin for t in t_tab])
    
   

'''Commande Constante
'''


def P_constante(t,P_constante):
    return P_constante
   
        
def P_constante_tab(t_tab,P):
    return np.array([P for t in t_tab])
    


'''Commande sinusoidale
'''


def P_sinusoidale(t,P_sin,T,phi,offset):
    
    res=offset+P_sin*np.sin(t/T*2*np.pi+phi)
    return res
   
        
def P_sinusoidale_tab(t_tab,P_sin,T,phi,offset):
    return np.array([ offset+P_sin*np.sin((t/T*2*np.pi)+phi) for t in t_tab])


def trace_plot(commande):
    sol = odeint(X_prime, [0, 0], t_tab, args = ( A, E0, I0, wEE, wEI, wII, wIE, alpha_E, alpha_I, w, phi,commande))
    
    if commande=="creneau":
        P=P_creneau_tab(t_tab,Pmin_creneau,Pmax_creneau,T_creneau)
    elif commande=="constante":
        P=P_constante_tab(t_tab,P_const)
    elif commande=="sinusoidale":
        P=P_sinusoidale_tab(t_tab,P_sin,T_sinusoidale,phi,offset)
    
    plt.figure()
    plt.plot()
    plt.plot(t_tab, sol[:, 0],'r')
    plt.plot(t_tab, P)
    
    plt.plot(t_tab, E0 + alpha_E * np.cos(w * t_tab))
    
    plt.savefig("Fig1")
   
    
    plt.figure()
    plt.plot()
    
    plt.plot(t_tab, P)
    
    plt.plot(t_tab, E0 + alpha_E * np.cos(w * t_tab))
    plt.plot(t_tab, sol[:, 1],'g')
    plt.savefig("Fig2")
    plt.show()


trace_plot("constante")
#trace_plot("sinusoidale")
#trace_plot("creneau")

