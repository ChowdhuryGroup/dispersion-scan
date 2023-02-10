# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:32:11 2021

@author: ChowdhuryGroup
"""
import csv
import matplotlib.pyplot as plt
backgroundname = "2021-06-29 dscan.txt"

with open(backgroundname) as bgfile:
    reader = csv.reader(bgfile,delimiter = '\t')

    plt.plot([float(i) for i in next(reader)],[float(i) for i in next(reader)])
    plt.show()