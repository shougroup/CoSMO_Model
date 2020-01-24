# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:49:43 2020

@author: ayuan
"""

import scipy.integrate
import scipy.interpolate
import numpy as np
import matplotlib.pyplot as plt

# L-A+ params
gL = 0.51; #max growth rate /hr
KL = 2.1; # Monod constant uM
nL = 3.2; # birth cooperativity
dL = 0.0024; # death rate /hr
cL = 5.4; # lysine consumption per birth fmole/cell
rA = 0.27; # hyp release rate fmole/cell/hr

#A-L+ params
gA = 0.44; # max growth rate /hr
KA = 1.3; # Monod constant uM
nA = 1.5; # birth cooperativity
dA = 0.015; # death rate /hr
cA = 3.1; # hyp consumption per birth fmole/cell
# lysine release rate interpolated in diffeq

#Sam's experiment 20150922

shour = [0.00, 5.83, 22.67, 29.50, 46.50, 53.50, 74.50]; # time (hrs)
#Live L-A+ density /mL, one experiment per row
sred = np.array([[1.75E6, 2.83E6, 4.95E6, 1.60E7, 3.03E8, 4.38E8, 1.69E9],
                 [1.75E6, 2.85E6, 4.88E6, 2.13E7, 2.85E8, 3.69E8, 1.19E9],
                 [1.75E6, 2.80E6, 4.72E6, 2.06E7, 2.90E8, 3.76E8, 1.55E9]]);
#Live A-L+ density /mL, one experiment per row
sgreen = np.array([[1.75E6, 1.81E6, 9.88E6, 8.89E6, 2.47E7, 2.78E8, 1.94E9],
                   [1.75E6, 1.80E6, 5.63E6, 5.95E6, 4.80E7, 3.54E8, 1.27E9],
                   [1.75E6, 1.86E6, 5.07E6, 5.69E6, 4.84E7, 3.61E8, 1.73E9]]);
#Total density /mL, one experiment per row
stotal = np.array([[3.50E6, 4.66E6, 1.50E7, 2.51E7, 3.28E8, 7.19E8, 3.64E9],
                   [3.50E6, 4.70E6, 9.86E6, 2.64E7, 3.39E8, 7.41E8, 3.30E9],
                   [3.50E6, 4.67E6, 1.06E7, 2.74E7, 3.34E8, 7.26E8, 2.47E9]]);

def cosmoFunc(y,t): # difeq for CoSMO growth and release
    NL, NA, SL, SA = y;
    Hyp = [0.00, 0.58, 0.65, 0.73, 0.80, 1.07]; # measured Hyp conc
    rLM = [0.52, 0.83, 0.78, 0.53, 0.38, 0.08]; # measured lys release rate
    fRel = scipy.interpolate.interp1d(Hyp,rLM);
    rLi = fRel(SA);
    
    if (rLi > 0):
        rL = rLi;
    else:
        rL = 0;

    dydt = [gL*SL**nL/(SL**nL + KL**nL)*NL - dL*NL,
            gA*SA**nA/(SA**nA + KA**nA)*NA - dA*NA,
            -cL*gL*SL**nL/(SL**nL + KL**nL)*NL + rL*NA,
            -cA*gA*SA**nA/(SA**nA + KA**nA)*NA + rA*NL];
            
    return dydt;

tstep = 0.5;
tspan = np.linspace(0,80,(80/tstep + 1));
y0 = [1.75E6/1000000, 1.75E6/1000000, 0, 0];
sol = scipy.integrate.odeint(cosmoFunc, y0, tspan);

plt.plot(tspan, sol[:,0], 'm', dashes=[3,2], label='NL (sim)')
plt.plot(tspan, sol[:,1], 'g', dashes=[3,2], label='NA (sim)')
plt.plot(tspan, sol[:,0] + sol[:,1], 'black', dashes=[3,2], label='Total (sim)')
plt.errorbar(shour, np.mean(sred,0)/1000000, 2*np.std(sred,0)/1000000, 
             color='m', label='NL (exp)')
plt.errorbar(shour, np.mean(sgreen,0)/1000000, 2*np.std(sgreen,0)/1000000,
             color='g', label='NA (exp)')
plt.errorbar(shour, np.mean(stotal,0)/1000000, 2*np.std(stotal,0)/1000000,
             color='black', label='Total (exp)')
plt.xlabel('Time (hours)')
plt.ylabel('Cell density (/ml)')
plt.legend(loc='best')
plt.yscale('log')
plt.show()