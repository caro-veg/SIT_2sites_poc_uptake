# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 11:02:22 2020

@author: Carolin
"""
import numpy as np
import subprocess


args = []
frange = np.linspace(0., 0.9, 10)
mutRates = np.array([2.45e-8, 4.9e-8, 1.23e-7, 2.45e-7, 2.45e-6, 2.45e-5, 7.95e-5])
for i in frange:
    for j in mutRates:
        arg = "6.02e-8 0.0054 0.0833 1.778 1.778 1 1 " + str(round(i, 2)) + " 2.45e-8 2.45e-8 " + str(j) + " " + str(j) + " 1478000 21780 0 220 0 1825 100"
        args.append(arg)

for i in frange:
    for j in mutRates:
        arg = "6.02e-8 0.0054 0.0833 1.778 1.778 1 1 " + str(round(i, 2)) + " 2.45e-8 2.45e-8 " + str(j) + " " + str(j) + " 1478000 21633 147 220 0 1825 100"
        args.append(arg)
        
for i in frange:
    for j in mutRates:
        arg = "6.02e-8 0.0054 0.0833 1.778 1.778 1 1 " + str(round(i, 2)) + " 2.45e-8 2.45e-8 " + str(j) + " " + str(j) + " 1478000 21120 660 220 0 1825 100"
        args.append(arg)

for i in frange:
    for j in mutRates:
        arg = "6.02e-8 0.0054 0.0833 1.778 1.778 1 1 " + str(round(i, 2)) + " 2.45e-8 2.45e-8 " + str(j) + " " + str(j) + " 1478000 18920 2860 220 0 1825 100"
        args.append(arg)

for i in frange:
    for j in mutRates:
        arg = "6.02e-8 0.0054 0.0833 1.778 1.778 1 1 " + str(round(i, 2)) + " 2.45e-8 2.45e-8 " + str(j) + " " + str(j) + " 1478000 13288 8492 220 0 1825 100"
        args.append(arg)



for i in range(0, len(args)):
    exe = ("C:\\Users\\cvegvari\\Documents\\Cpp\\SIT_2Sites\\TreatmentScenarios\\SIT_2sites_poc_uptake\\bin\\Release\\SIT_2sites_poc_uptake.exe" + " " + args[i]).split()
    popen = subprocess.Popen(exe)
    popen.wait()
    print("finished!" + " " + str(i))
    


