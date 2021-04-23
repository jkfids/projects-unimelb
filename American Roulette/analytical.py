# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 17:52:53 2021

@author: Fidel
"""

from numpy import log

f = 0.05
q = 9/19
tau = 60

x_w = 2000
x_m = 10
x_0 = 1000

def alpha(f, q, tau):
    return -2*(1-2*q)/f

def beta(f, q, tau):
    return -2*tau/(f**2)
    
alph = alpha(f, q, tau)
bet = beta(f, q, tau)

def c_1(x_w, x_m, f=f, q=q, tau=tau):
    a = 1-alpha(f, q, tau)
    b = beta(f, q, tau)
    return b*(log(x_w)-log(x_m))/(a*(x_w**a-x_m**a))

c1 = c_1(x_w, x_m)

def c_0(x_w, c1, f=f, q=q, tau=tau):
    a = 1-alpha(f, q, tau)
    b = beta(f, q, tau)
    return b*log(x_w)/a - c1*x_w**a

c0 = c_0(x_w, c1)

print(c0, c1)

def calc_T(x_0, c0, c1, f=f, q=q, tau=tau):
    a = 1-alpha(f, q, tau)
    b = beta(f, q, tau)
    return c0 + c1*x_0**a - b*log(x_0)/a
    
T = calc_T(x_0, c0, c1)
print(T)