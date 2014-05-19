#!/usr/bin/python

import matplotlib.pyplot as plt
from numpy.random import *

import etframes


data = uniform(0, 99, 10) + 1.0

etframes.bar_chart(data, [25, 50, 75])
plt.xticks(range(len(data)))

plt.show()

