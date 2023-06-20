#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 15:37:05 2022

@author: laidlaw
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

k1 = 1
k2 = 1
q1 = -1
q2 = 1
q3 = 1
q4 = -1
l1 = 2
l2 = 2

m1 = 1
m2 = 1
m3 = 2
m4 = 2


def Pot(x1, x2, x3, x4):
    temp = 0.5*k1*(x2-x1-l1)**2
    temp += 0.5*k2*(x4-x3-l2)**2
    temp += q1*q2/(x2-x1)
    temp += q1*q3/(x3-x1)
    temp += q1*q4/(x4-x1)
    temp += q2*q3/(x3-x2)
    temp += q2*q4/(x4-x2)
    temp += q3*q4/(x4-x3)
    return temp

dt = 0.001
dx = 0.000001
length = 50000

X1 = np.zeros(length)
X2 = np.zeros(length)
X3 = np.zeros(length)
X4 = np.zeros(length)
V1 = np.zeros(length)
V2 = np.zeros(length)
V3 = np.zeros(length)
V4 = np.zeros(length)
time = np.zeros(length)

X1[0] = -21
X2[0] = -20
X3[0] = 20
X4[0] = 21

V1[0] = 5
V2[0] = 5
V3[0] = -5
V4[0] = -5

for i in range(len(X1)-1):
    X1[i+1] = X1[i] + V1[i]*dt
    X2[i+1] = X2[i] + V2[i]*dt
    X3[i+1] = X3[i] + V3[i]*dt
    X4[i+1] = X4[i] + V4[i]*dt
    a1 = -1/m1*(Pot(X1[i]+dx,X2[i],X3[i],X4[i])-Pot(X1[i],X2[i],X3[i],X4[i]))/dx
    a2 = -1/m2*(Pot(X1[i],X2[i]+dx,X3[i],X4[i])-Pot(X1[i],X2[i],X3[i],X4[i]))/dx
    a3 = -1/m3*(Pot(X1[i],X2[i],X3[i]+dx,X4[i])-Pot(X1[i],X2[i],X3[i],X4[i]))/dx
    a4 = -1/m4*(Pot(X1[i],X2[i],X3[i],X4[i]+dx)-Pot(X1[i],X2[i],X3[i],X4[i]))/dx
    V1[i+1] = V1[i] + a1*dt
    V2[i+1] = V2[i] + a2*dt
    V3[i+1] = V3[i] + a3*dt
    V4[i+1] = V4[i] + a4*dt
    time[i+1] = time[i] + dt

plt.plot(time,X1)
plt.plot(time,X2)
plt.plot(time,X3)
plt.plot(time,X4)
#plt.axis([0,8,-10,10])
plt.show()

plt.plot(time, V1)
plt.plot(time, V2)
plt.plot(time, V3)
plt.plot(time, V4)
plt.show()