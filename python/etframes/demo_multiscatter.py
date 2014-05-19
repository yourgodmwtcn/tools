#!/usr/bin/python

import matplotlib.pyplot as plt
from numpy.random import *

import etframes

data = normal(size=(100,4))
data[:,2] += data[:,1] # lets add a linear dependency

etframes.multi_scatter(data, ('V1', 'V2', 'V3', 'V4'))

plt.show()

