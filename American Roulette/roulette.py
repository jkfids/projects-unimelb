# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 20:34:31 2021

@author: Fidel
"""

import numpy as np
from random import random
from matplotlib import pyplot as plt
from time import time
from numpy.polynomial import Polynomial

f = 0.05
q = 9/19
tau = 60

x_0 = 1000
x_w = 2000
x_m = 10

class Roulette:
    def __init__(self, x_0=x_0, x_w=x_w, x_m=x_m, f=f, q=q, tau=tau):
        self.f = f
        self.q = q
        self.tau = tau
        self.x = x_0
        self.x_w = x_w
        self.x_m = x_m
        self.tau = tau
        self.t = 0
        
    def play(self):
        rand = random()
        if rand < self.q:
            self.x += self.f*self.x
        elif rand > self.q:
            self.x -= self.f*self.x
        self.t += self.tau
            
    def start(self):
        while (self.x > self.x_m) & (self.x < self.x_w):
            self.play()

def plot_t_pdf(N, x_0=x_0, x_w=x_0, x_m=x_m, f=f, q=q, tau=tau):
    start = time()
    fig, ax = plt.subplots(dpi=144)
    
    t_history = np.zeros(N, dtype=np.int32)
    for i in range(N):
        andrew = Roulette()
        andrew.start()
        t_history[i] = andrew.t
        
    t_mean = np.mean(t_history)
    bin_heights, bin_borders, _ = ax.hist(t_history, density=True, bins='auto')
    bin_centers = bin_borders[:-1] + np.diff(bin_borders)/2
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Probability')
    end = time()
    fig.savefig('playtime2_pdf')
    print(f'Time elapsed: {round(end-start, 2)}s')  
    return t_mean, t_history, bin_centers, bin_heights

t_mean, t_history, bin_centers, bin_heights = plot_t_pdf(100000)
print(f'Average playing time: {round(t_mean)}s')
