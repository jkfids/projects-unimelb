#!/usr/bin/env python3
"""Sandpile Lab
============

This is the main file from which we run all the various routines in the
sandpile lab.

"""
from pathlib import Path

import numpy as np
import scipy as sp
import scipy.stats
from matplotlib import pyplot

from sandpile import SandPile

"""Linear function used for fitting a line of best fit"""
def lin_func(x, m, c):
    return m*x + c

def massovertime_centre_plot(output_dir):
    """Plots mass on grid over time for a 3 by 3 matrix when grains of 
       are dropped in the centre of the grid."""
       
    massvtime_centre = SandPile(3,3)
    for t in range(100):
        massvtime_centre.drop_sand(1,[1,1]);
    
    x = [i for i in range(101)];
    y = massvtime_centre.mass_history;

    fig, ax = pyplot.subplots()
    ax.plot(x, y)
    ax.set_title("Mass on Grid over Time for Centre Dropped Grains")
    ax.set_xlabel("Time")
    ax.set_ylabel("Mass on Grid")

    fig.savefig("output/massovertime_centre_plot.jpg")
    pyplot.close(fig)

def massovertime_random_plot(output_dir):
    """Plots mass on grid over time for a 3 by 3 matrix when grains of 
       are dropped randomly."""
       
    massvtime_random = SandPile(3,3);
    for t in range(100):
        massvtime_random.drop_sand();
        
    x = [i for i in range(101)];
    y = massvtime_random.mass_history;

    fig, ax = pyplot.subplots()
    ax.plot(x, y)
    ax.set_title("Mass on Grid over Time for Randomly Dropped Grains")
    ax.set_xlabel("Time")
    ax.set_ylabel("Mass on Grid")

    fig.savefig("output/massovertime_random_plot.jpg")
    pyplot.close(fig)
    


def massovertime_scale_plot(output_dir):
    """Compares plots of mass on grid over time for a 3x3 matrix, 9x9 matrix
       15x15 matrix for randomly dropped grains."""
       
    massvtime_9 = SandPile(9,9)
    massvtime_15 = SandPile(15,15)
    massvtime_21 = SandPile(21,21)
    for t in range(1250-1):
        massvtime_21.drop_sand();
        massvtime_9.drop_sand();
        massvtime_15.drop_sand();
    
    x = [i for i in range(1250)];
    y1 = (massvtime_9.mass_history);
    y2 = (massvtime_15.mass_history);
    y3 = (massvtime_21.mass_history);

    fig, ax = pyplot.subplots()
    ax.plot(x, y1, label='9 by 9 Grid')
    ax.plot(x, y2, label='15 by 15 Grid')
    ax.plot(x, y3, label='21 by 21 Grid')
    #ax.set_title("Mass on Grid over Time across Grid Sizes - 1250 iterations")
    ax.set_xlabel("Time")
    ax.set_ylabel("Mass on Grid")
    ax.legend()

    fig.savefig("output/massovertime_scale.jpg")
    pyplot.close(fig)
    
def densityovertime_centre_plot(output_dir):
    """Plots the mass density on the grid over time for a 9 by 9 matrix for
    centre dropped grains over 1000 iterations, and calculates the average
    grain density over a certain period of time"""
    
    densitygrid = SandPile(9,9)
    for t in range(1000-1):
        densitygrid.drop_sand(1,[4,4])
    
    x = [i for i in range(1000)];
    y = np.divide(densitygrid.mass_history,81)
    
    fig, ax = pyplot.subplots()
    ax.plot(x, y)
    ax.set_title("Density on Grid over Time (Centre) - 1000 iterations, 9x9 grid")
    ax.set_xlabel("Time")
    ax.set_ylabel("Mass Density of Grid")
    fig.savefig("output/densityovertime_centre_plot.jpg")
    pyplot.close(fig)
    
    y_average = sum(y[200:1000])/800;
    print('Average grain density (centre) of 9x9 grid over 200<t<1000:', y_average)
    
def densityovertime_random_plot(output_dir):
    """Plots the mass density on the grid over time for a 9 by 9 matrix for
    randomly dropped grains over 1000 iterations, and calculates the average
    grain density over a certain period of time"""
    
    density_random = SandPile(9,9)
    for t in range(1000-1):
        density_random.drop_sand()
    
    x = [i for i in range(1000)];
    y = np.divide(density_random.mass_history,81)
    
    fig, ax = pyplot.subplots()
    ax.plot(x, y)
    ax.set_title("Density on Grid over Time (Random) - 1000 iterations, 9x9 grid")
    ax.set_xlabel("Time")
    ax.set_ylabel("Mass Density of Grid")
    fig.savefig("output/densityovertime_random_plot.jpg")
    pyplot.close(fig)
    
    y_average = sum(y[200:1000])/800;
    print('Average grain density (random) of 9x9 grid over 200<t<1000:', y_average)
    
def topplesovertime_plot(output_dir):
    """Plots the number of topples per avalanche over time for a 10x10 grid
    across 600 iterations."""
    
    topplesovertime = SandPile(10,10)
    for t in range(1000):
        topplesovertime.drop_sand()
    
    x = [i for i in range(1000)];
    y = topplesovertime.topples_history;
    
    fig, ax = pyplot.subplots()
    ax.plot(x, y)
    ax.set_title("Number of Topples over Time - 500 iterations, 10x10 grid")
    ax.set_xlabel("Time")
    ax.set_ylabel("Number of Topples")
    
    fig.savefig("output/topplesovertime_plot.jpg")
    pyplot.close(fig)

"""
10x10 sandpile simulated over 20000 iterations in order to graph the
frequency of various avalanche properties
"""

avalancheproperties = SandPile(10,10)
for t in range (20000-1):
    avalancheproperties.drop_sand()    

def freqvtopples_plot(output_dir):
    """Plots the frequency of number of topples per avalanche for 
    20,000 iterations in a 10x10 grid"""

    history = np.unique(avalancheproperties.topples_history,return_counts=True)
    x = history[0];
    x = np.delete(x,0)
    y = history[1];
    y = np.delete(y,0)
    ln_y = np.log(y)

    fig, ax = pyplot.subplots()
    ax.scatter(x, ln_y)
    ax.set_title("Frequency of No. of Topples - 20000 iterations, 10x10 grid")
    ax.set_xlabel("Number of Topples")
    ax.set_ylabel("ln(Frequency)")

    fig.savefig("output/freqvtopples_plot.jpg")
    pyplot.close(fig)
    
def freqvtopples_loglog_plot(output_dir):
    """Plots the frequency of number of topples per avalanche for 
    20,000 iterations in a 10x10 grid, on a log-log scatter plot. Also applies
    a line of best fit over the linear region."""

    history = np.unique(avalancheproperties.topples_history,return_counts=True)
    x = history[0];
    x = np.delete(x,0)
    ln_x = np.log(x)
    y = history[1];
    y = np.delete(y,0)
    ln_y = np.log(y)
    
    # Fits the data to a line of best fit using sp.optimize.curve_fit
    # and the lin_func we previously defined.
    ln_x_range = ln_x[ln_x<3.5]
    ln_y_range = ln_y[0:len(ln_x_range)]
    param = sp.optimize.curve_fit(lin_func, ln_x_range, ln_y_range)
    [m, c] = param[0]
    y_fit = (m*ln_x_range+c)
    # The two lines below extend the line fit
    ln_x_range = np.append(ln_x_range,5)
    y_fit = np.append(y_fit,m*5+c)
    
    fig, ax = pyplot.subplots()
    ax.scatter(ln_x, ln_y)
    string = "ln(y) = {}ln(x) + {}"
    string = string.format(round(m,5),round(c,5))
    ax.plot(ln_x_range, y_fit, label=string, c='darkorange', ls = '--')
    
    #ax.set_title("Frequency of No. of Topples - 20000 iterations, 10x10 grid")
    ax.set_xlabel("ln(Number of Topples)")
    ax.set_ylabel("ln(Frequency)")
    ax.legend()

    fig.savefig("output/freqvtopples_loglog_plot.jpg")
    pyplot.close(fig)
    
def freqvarea_plot(output_dir):
    """Compares the area affected per avalanche for 20,000 iterations in a 
    10x10 grid"""

    history = np.unique(avalancheproperties.area_toppled_history,return_counts=True)
    x = history[0];
    x = np.delete(x,0)
    y = history[1];
    y = np.delete(y,0)
    ln_y = np.log(y)

    fig, ax = pyplot.subplots()
    ax.scatter(x, ln_y)
    ax.set_title("Frequency of Avalanche Area - 20000 iterations, 10x10 grid")
    ax.set_xlabel("Area of Avalanche")
    ax.set_ylabel("ln(Frequency)")

    fig.savefig("output/freqvarea_plot.jpg")
    pyplot.close(fig)
    
def freqvloss_plot(output_dir):
    """Compares the area affected per avalanche for 20,000 iterations in a 
    10x10 grid"""

    history = np.unique(avalancheproperties.grain_loss_history,return_counts=True)
    x = history[0];
    x = np.delete(x,0)
    y = history[1];
    y = np.delete(y,0)
    ln_y = np.log(y)

    fig, ax = pyplot.subplots()
    ax.scatter(x, ln_y)
    ax.set_title("Frequency of Grain Loss per Avalanche - 20000 iterations, 10x10 grid")
    ax.set_xlabel("Grain Loss per Avalanche")
    ax.set_ylabel("ln(Frequency)")

    fig.savefig("output/freqvloss_plot.jpg")
    pyplot.close(fig)
    
def freqvMlength_plot(output_dir):
    """Compares the frequency of avalanche Manhattan lengths for 20,000 
    iterations in a 10x10 grid"""

    history = np.unique(avalancheproperties.M_length_history,return_counts=True)
    x = history[0];
    x = np.delete(x,0)
    ln_x = np.log(x)
    y = history[1];
    y = np.delete(y,0)
    ln_y = np.log(y)

    fig, ax = pyplot.subplots()
    ax.scatter(x, ln_y)
    ax.set_title("Frequency of Avalanche Lengths - 20000 iterations, 10x10 grid")
    ax.set_xlabel("Avalanche Manhattan Length")
    ax.set_ylabel("ln(Frequency)")

    fig.savefig("output/freqvMlength_plot.jpg")
    pyplot.close(fig)
    
def freqvElength_plot(output_dir):
    """Compares the frequency of avalanche Euclidean lengths for 20,000 
    iterations in a 10x10 grid"""

    history = np.unique(avalancheproperties.E_length_history,return_counts=True)
    x = history[0];
    x = np.delete(x,0)
    y = history[1];
    y = np.delete(y,0)
    ln_y = np.log(y)

    fig, ax = pyplot.subplots()
    ax.scatter(x, ln_y)
    ax.set_title("Frequency of Avalanche Lengths - 20000 iterations, 10x10 grid")
    ax.set_xlabel("Avalanche Euclidean Length")
    ax.set_ylabel("ln(Frequency)")

    fig.savefig("output/freqvElength_plot.jpg")
    pyplot.close(fig)
    
def areavlength_plot(output_dir):
    """Plots the correlation between avalanche area and length
    on a log-log scatter plot for a 10x10 grid across 20,000 iterations."""

    x = np.array(avalancheproperties.M_length_history);
    x = x[x>0]
    ln_x = np.log(x)
    y = np.array(avalancheproperties.area_toppled_history);
    y = y[y>0]
    ln_y = np.log(y)

    fig, ax = pyplot.subplots()
    ax.scatter(ln_x, ln_y)
    ax.set_title("Avalanche Area vs Length - 20000 iterations, 10x10 grid")
    ax.set_xlabel("ln(Avalanche Manhattan Length)")
    ax.set_ylabel("ln(Avalanche Area)")
    
    fig.savefig("output/areavlength_plot.jpg")
    pyplot.close(fig)
    
def areavtopples_plot(output_dir):
    """Plots the correlation between avalanche area and number of topples
    on a log-log scatter plot for a 10x10 grid across 20,000 iterations."""

    x = np.array(avalancheproperties.topples_history);
    y = np.array(avalancheproperties.area_toppled_history);

    fig, ax = pyplot.subplots()
    ax.scatter(x, y)
    ax.set_title("Avalanche Area vs Number of Topples - 20000 iterations, 10x10 grid")
    ax.set_xlabel("Number of Topples")
    ax.set_ylabel("Avalanche Area")
    
    fig.savefig("output/areavtopples_plot.jpg")
    pyplot.close(fig)
    

    
"""
Square sandpile grids of increasing size iterated over 100000 time steps. 
These are commented out since they require a long amount of time to simulate.
If you have 5 minutes to spare you can uncomment up to grid size 45. If you
have half an hour and a decent computer you can uncomment everything.
"""

data = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
fig, ax = pyplot.subplots()
ax.imshow(data)
ax.set_title("Grid plot")
fig.savefig("output/example_grid.jpg")
pyplot.close(fig)


"""square9 = SandPile(9,9)
square9.drop_sand(100000, [4,4])
fig, ax = pyplot.subplots()
ax.axis('off')
ax.imshow(square9.grid)
ax.set_title("9 x 9 Grid")
fig.savefig("output/square9_grid.jpg")
pyplot.close(fig)

square15 = SandPile(15,15)
square15.drop_sand(100000, [7,7])
fig, ax = pyplot.subplots()
ax.imshow(square15.grid)
ax.set_title("15x15 Grid Plot - 100000 iterations")
fig.savefig("output/square15_grid.jpg")
pyplot.close(fig)

square27 = SandPile(27,27)
square27.drop_sand(100000, [13,13])
fig, ax = pyplot.subplots()
ax.imshow(square27.grid)
ax.set_title("27x27 Grid Plot - 100000 iterations")
fig.savefig("output/square27_grid.jpg")
pyplot.close(fig)

square45 = SandPile(45,45)
square45.drop_sand(100000, [22,22])
fig, ax = pyplot.subplots()
ax.imshow(square45.grid)
ax.set_title("45x45 Grid Plot - 100000 iterations")
fig.savefig("output/square45_grid.jpg")
pyplot.close(fig)

square69 = SandPile(69,69)
square69.drop_sand(100000, [34,34])
fig, ax = pyplot.subplots()
ax.axis('off')
ax.imshow(square69.grid)
ax.set_title("69 x 69 Grid")
fig.savefig("output/square69_grid.jpg")
pyplot.close(fig)

square99 = SandPile(99,99)
square99.drop_sand(100000, [49,49])
fig, ax = pyplot.subplots()
ax.axis('off')
ax.imshow(square99.grid)
ax.set_title("99 x 99 Grid")
fig.savefig("output/square99_grid.jpg")
pyplot.close(fig)

square135 = SandPile(135,135)
square135.drop_sand(100000, [67,67])
fig, ax = pyplot.subplots()
ax.axis('off')
ax.imshow(square135.grid)
ax.set_title("135 x 135 Grid")
fig.savefig("output/square135_grid.jpg")
pyplot.close(fig)"""


def main():
    # Make sure that the output/ directory exists, or create it otherwise.
    output_dir = Path.cwd() / "output"
    if not output_dir.is_dir():
        output_dir.mkdir()

    massovertime_centre_plot(output_dir)
    massovertime_random_plot(output_dir)
    massovertime_scale_plot(output_dir)
    densityovertime_centre_plot(output_dir)
    densityovertime_random_plot(output_dir)
    topplesovertime_plot(output_dir)
    
    freqvtopples_plot(output_dir)
    freqvtopples_loglog_plot(output_dir)
    freqvarea_plot(output_dir)
    freqvloss_plot(output_dir)
    freqvMlength_plot(output_dir)
    freqvElength_plot(output_dir)
    areavlength_plot(output_dir)
    areavtopples_plot(output_dir)
    
    
    

if __name__ == "__main__":
    main()
