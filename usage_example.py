"""
Usage example for rga_data_handling code package

@author: Georg Schlisio
@contact: georg.schlisio@ipp.mpg.de
@date: 2019-11-27
"""

import numpy as np
import matplotlib.pyplot as pp
import random

import rga_data_handling.Class as Class
from rga_data_handling.RGA_class_plotter import plot_data, show_results
import rga_data_handling.fit_many_times as fmt

# synthetic test data
columns = {
    "time": np.linspace(0, 100, 100),
    1: np.array([random.random() for i in range(100)]) * 1e-6,
    2: np.array([random.random() for i in range(100)]) * 1e-7,
    4: np.array([random.random() for i in range(100)]) * 1e-8,
}
# annotation of columns
tags = {"time_col": "time", "title": "testdata"}

# create object from data and annotations
tr = Class.Trace(columns, tags)

# Plotting the given data
# Optionally give subset of signals as second argument, i.e. masses 1 and 3:
# plot_data(tr, [1, 3])
# To save the plot to file:
# plot_data(tr, save_name='plot_image.png')
plot_data(tr)

# non-zero background before the pulse? (e.g. from -15 s to -5 s)
# tr.subtract_background([-15, -5])

# Offset correction
# Correcting for offsets in channels where no signal is expected and signal is
# considered noise (e.g. from magnetic field, pressures jumps etc)
# tr.correct_offset(23)
# or multiple (taking average signal as offset)
# tr.correct_offset([23, 37, 43])

## Define relevant molecules
# Molecules are split into two groups: containing hydrogen and not containing hydrogen
# both groups are defined with a speaking name, a unique fitting parameter name and
# the molecule name according to the library (predefined names and values can be found in
# RGA_calibration_reader.py).
#
# For the hydrogen-containing species, the H/(H+D) ratio can also be given as either fixed
# value between 0 and 1 or as unique fitting parameter name.
# The H/(H+D) ratio can also be restricted as shown below.
#
# ATTENTION: Use only letters and no whitespaces for unique fitting parameter names, since
# they will be used in formulas which need to be evaluated!
# A re-use of fitting parameters will couple the variables!

# hydrogen-containing molecules
H_molecules = {
    "water": ["w", "water", "rw"],  # water molecule with free H/(H+D) ratio
    "methane": ["m", "methane", 0.8],  # methane molecule with fixed H/(H+D) ratio
    "ammonia": [
        "a",
        "ammonia",
        ["ra", 0, 0.1],
    ],  # ammonia molecule with restricted H/(H+D) ratio
}

# hydrogen-free molecules
non_H_molecules = {
    "n2": ["n", "N2"],
    "o2": ["o", "O2"],
}

# conditionally ignoring masses
# ignore complete signal by giving its name or ignore signal if corresponding signal is above threshold
# e.g. at high N2 pressures, the signal at 14amu is dominated by N2. Because N2 is sufficienlty well identified by 28amu,
# looking at 14amu too just causes unnecessary error in the fit
disregard = [
    16,  # completely leave out 16amu
    [14, 28, 1e-7],  # ignore 14amu if 28amu exceeds 1e-7
]

# restrict evaluation to given interval and timestep
tstart = 30
tstop = 50
tstep = 1

# deconvolution of data - fitting
tr.deconvolute(H_molecules, non_H_molecules, disregard, tstart, tstop, tstep)

# show deconvolution results, optionally restricted to some molecules
show_results(tr)
show_results(tr, ["water", "methane", "ammonia"])

# what is the effect of poorly defined cracking patterns?
# perform the fit many times, each time with perturbed cracking patterns, then look at the average value

# To assess influence of cracking pattern uncertainties, make perturbations and fit with different combinations
# Note: the fourth argument in fmt is a relict from early versions that has not yet been moved to a proper place
perturb = 0.1  # 10 % error of the cracking pattern
times = 25  # number of different perturbations to be tested
traces = fmt.fit_many_times(
    tr,
    times,
    perturb,
    [[], ""],
    H_molecules,
    non_H_molecules,
    disregard,
    tstart,
    tstop,
    tstep,
)

# average over result traces
average_trace = fmt.averages(traces)

# plot results with error bars
fmt.errorbar_results(average_trace)

# show plots
pp.show()
