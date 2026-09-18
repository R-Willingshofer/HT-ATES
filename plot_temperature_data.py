# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 09:54:56 2026

@author: robin
"""
#import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df_DARTS = pd.read_excel('Preliminary_synthetic_data_DARTS_r0.xlsx')
df_Seawat = pd.read_csv('tempdataR0.csv')
df_fimbull = pd.read_csv('daily_aquifer_averaged_temperatures.csv')
df_MF = pd.read_csv('mf6_well_temperature.csv')

time_Seawat = df_Seawat['Unnamed: 0'].values
temp_Seawat = df_Seawat['temp'].values

time_Darts = df_DARTS['time (days)'].values
temp_Darts = df_DARTS['Temperature (C)'].values

time_MF = df_MF['time_d'].values
temp_MF = df_MF['T_r000'].values

time_FB = df_fimbull['time_d'].values
temp_FB = df_fimbull['T_r000'].values

nyears = 5
years = np.arange(0, nyears, 1)
days_year = 360

start_charge = (years) * days_year
stop_charge = (years + (90/360)) * days_year

start_discharge = (years + (180/360)) * days_year
stop_discharge = (years + (270/360)) * days_year

plt.figure(figsize = (7, 4))

for i in range(nyears):
    plt.axvspan(start_charge[i], stop_charge[i], facecolor='b', alpha=0.1, zorder = 0)
    plt.axvspan(start_discharge[i], stop_discharge[i], facecolor='r', alpha=0.1, zorder = 0)


plt.plot(time_Darts, temp_Darts, label = 'DARTS')
plt.plot(time_Seawat, temp_Seawat, label = 'Seawat')
plt.plot(time_MF, temp_MF, label = 'Modflow')
plt.plot(time_FB, temp_FB, label = 'fimbull')
plt.ylabel('Temperature (C)')
plt.xlabel('Time (d)')
plt.xlim((0, 1800))
plt.grid()
plt.legend()
plt.show()