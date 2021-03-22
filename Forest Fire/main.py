# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 16:24:11 2021

@author: Fidel
"""

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.animation import FuncAnimation
from time import time

from forestfire import ForestFire

def animate_forest(forest, interval=100, frames=100, name='forestfire.gif'):
    # Animate a forest fire for a given number of frames (i.e. timesteps)
    start = time()
    cmap = colors.ListedColormap(['red', 'black', 'green'])
    bounds = [-1, -0.5, 0.5, 1]
    norm = colors.BoundaryNorm(bounds,  cmap.N)
    
    fig, ax = plt.subplots()
    ax.axis('off')
    
    fig = plt.figure(frameon=False)
    fig.set_size_inches(10,10)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    
    def init_frame():
        ax.imshow(forest.grid, cmap=cmap, norm=norm, aspect='auto')
    
    def animate(i):
        ax.imshow(forest.grid, cmap=cmap, norm=norm, aspect='auto')
        forest.step()
        print(i)
    
    anim = FuncAnimation(fig, animate, init_func=init_frame, interval=interval, frames=frames)
    anim.save('animations/' + name)
    end = time()
    print(f'Time elapsed: {round((end - start), 2)} seconds')
    
def plot_svt(forest, t_max):
    # Plot fire size vs t
    fig, ax = plt.subplots()
    ax.set_xlabel('Time')
    ax.set_ylabel('Fire Size')

    props = dict(boxstyle='square', facecolor='white')
    textbox = (
                f'L = {forest.height}\n'
                f'p = {forest.p}\n'
                f'f = {forest.f}'
    )
    ax.text(0.025, 0.965, textbox, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
    
    forest.step(t_max)
    y = forest.s_history
    x = range(len(y))
    
    ax.plot(x, y)

def pd_firesize(forest, t, N):
    # Constructs a histogram of probability vs. fire size after certain time
    start = time()
    L = forest.width
    f = forest.f
    p = forest.p
    
    firesizes = []
    for i in range(N):
        forest = ForestFire([L,L], f, p, spark=True)
        forest.step(t)
        if forest.count(-1) != 0:
            firesizes.append(forest.count(-1))
        #print(i)
    plt.hist(firesizes, density=True, bins=35)
    plt.ylabel('Probability')
    plt.xlabel('Fire Size')
    plt.title('Fire Size Probability Distribution')
    end = time()
    print(f'Time elapsed: {round((end - start), 2)} seconds')
     
#%%
L = 100
high_grow_no_lightning = ForestFire([L,L], 0, 0.5)
pd_firesize(high_grow_no_lightning, 150, 100)
    
#%%
L = 200
forest_100=ForestFire([L,L], 0, 0.5, spark=True)
plot_svt(forest_100, 200)

#%%
L = 720
forest = ForestFire([L,L], 0.0001, 0.01)
animate_forest(forest, interval=100, frames=200)