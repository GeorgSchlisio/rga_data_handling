# RGA data handling and analysis code

## Purpose
Mass spectrometry in fusion environments involves multiple species of hydrogen.
This code allows to treat (partially) deuterated species as commonly found in current fusion experiments working with H2 and D2.
It provides a set of functions for handling of RGA recordings, viewing and analysis.

For a quick start read `usage_example.py`.

Further reading: https://doi.org/10.1016/j.fusengdes.2017.05.037

## Installation
Clone repository

`git clone `

Checkout `refactor` branch

`git checkout refactor`

To install within Python setup

`python setup.py install`

## Dependencies
* `sympy`
* `matplotlib.pyplot`
* `scipy`
* `os.path`
* `time`

## Code structure
Class.py: RGA data objects, related functions
molecules2.py: Functions for handling of cracking patterns
RGA_fitting.py: Functions for fitting cracking patterns on RGA recordings
RGA_class_plotter.py: Functions for viewing data and results

## Authors
Original Author:
Aleksander Drenik:
IPP, Garching;
IJS, Ljubljana

Maintainer of this fork:
Georg Schlisio:
IPP, Greifswald

