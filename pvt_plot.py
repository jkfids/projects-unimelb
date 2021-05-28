# -*- coding: utf-8 -*-
"""
Created on Fri May 21 12:56:02 2021

@author: Fidel
"""

import numpy as np
from numpy import pi, exp
import matplotlib.pyplot as plt

def myfn(t, w=1):
    return (t*t/w)*exp(-w*w*t*t/2)

X = np.linspace(0, 12, 1001)
Y1 = myfn(X, 1)
Y2 = myfn(X, 0.8)
Y3 = myfn(X, 0.6)
Y4 = myfn(X, 0.4)
   
plt.plot(X, Y1, label='\u03C9 = 1.0')
plt.plot(X, Y2, label='\u03C9 = 0.8')
plt.plot(X, Y3, label='\u03C9 = 0.6')
plt.plot(X, Y4, label='\u03C9 = 0.4')
plt.xlabel('\u03C4') 
plt.ylabel('Probability')
plt.legend()


ax = plt.gca()
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

plt.grid(True)
plt.savefig('pvt_plot.png', dpi=144)