# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 14:19:50 2024

@author: Titus
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import math
#%% Import constants and intial conditions

t_half = pd.read_csv('Input/Half_life.csv') #half-life in years of isotopes of intrest
g_not = pd.read_csv('Input/Initial_Isotope_inventory.csv', header=1) # intial inventory of isotopes of intrest

N_not = g_not['initial inventory'] / t_half['g/mol'] # convert grams to moles
lamda = np.divide( np.log(2), t_half['half-life (yr)'] ) # calcualtes the exponential decay constant lamda = ln(2)/halt-life
assume = pd.read_csv('Input/groundwater_and_dissolution.csv', header=2) # import other assumed values

#%% Import ultrasonic monitor data
UT_df = pd.read_csv('Input/UT_data.csv',header=2)
#%% calculate canister thickness and corrosion rate

V = 5949 -0.9723 * UT_df.loc[0,'T'] #logitudnal speed of sound in carbon steel T is in °C

D = UT_df.loc[0,'t_eco1'] * V / 2

CorRate = (UT_df.loc[0,'t_eco2'] - UT_df.loc[0,'t_eco1'])* V / (2 * UT_df.loc[0,'delta_t_eco'])  

t_fail = (D - assume.loc[0,'D_min']) / CorRate + UT_df.loc[0,'t_past']

#%% Decay chain 1 is alpha decay of Am-243
t=100 #time in years

def Am243decay(N_not, t): # N_not is a series where N_not[3] is the number of moles of Am-243 at t=0
    Nt_Am243 = N_not[3] * np.e ** (-lamda[3] * t) #Bateman equation solution for single decay
    return Nt_Am243

#%% 4n+1 decay chain Cm-245 -> Pu-241 -> Am-241

def Am241chain(N_not, t): # N_not had the intial moles (in this order) of Cm-245, Pu-241, Am-241
    Nt_Cm245 = N_not[0] * np.e ** (-lamda[0] * t)
    b = N_not[0] * lamda[0] / (lamda[1] - lamda[0])
    Nt_Pu241 = b * np.e ** (-lamda[0] *t) - (b * b * np.e ** (-lamda[1] *t)) + N_not[1] * np.e ** (-lamda[1] *t) 
    Nt_Am241 = lamda[1] *b / (lamda[2] - lamda[0]) *(np.e ** (-lamda[0] *t) - np.e ** (-lamda[2] *t)) - \
        lamda[1] * b / (lamda[2] - lamda[1]) *(np.e ** (-lamda[1] *t) - np.e ** (-lamda[2] *t)) + \
        lamda[1] * N_not[1] / (lamda[2] - lamda[1]) *(np.e ** (-lamda[1] *t) - np.e ** (-lamda[2] *t)) + \
        N_not[2] * np.e ** (-lamda[2] * t)
    return Nt_Cm245, Nt_Pu241, Nt_Am241

#%% Calculate decay pofiles

t_steps = np.linspace(0, 1e6, num=500) # set times to record isotope quantites at

#Make data frame to hold the concentration profile data
n_timesteps = range(len(t_steps))
N_t = pd.DataFrame(data=t_steps, columns=["Time (yr)"])
N_t['Am-243']= np.zeros(len(t_steps))
N_t['Cm-245']= np.zeros(len(t_steps))
N_t['Pu-241']= np.zeros(len(t_steps))
N_t['Am-241']= np.zeros(len(t_steps))

#Calculate and store to df the isotope concentration profiles
for x in n_timesteps:
   
    N_t.iat[x,1] = Am243decay(N_not, t_steps[x])
    N_t.iat[x,2], N_t.iat[x,3], N_t.iat[x,4] = Am241chain(N_not, t_steps[x])

# Calculate and store to df the isotope concentration at time of canister failer    
N_fail = pd.DataFrame(g_not['Isotope'])
N_fail['moles'] = np.zeros(4)
N_fail['Time (yr)']  =np.full(4, t_fail)  

N_fail.loc[3,'moles'] = Am243decay(N_not, t_fail)
N_fail.loc[0,'moles'], N_fail.loc[1,'moles'], N_fail.loc[2,'moles'] = Am241chain(N_not, t_fail)
#%% plot decay profile

plt.style.use('Input/publish3.mplstyle')
fig, ax = plt.subplots(figsize=(5,5)) #size is in inches
ax.plot(N_t["Time (yr)"], N_t['Am-243'], label='Am-243')
ax.plot(N_t["Time (yr)"], N_t['Cm-245'], label='Cm-245')
ax.plot(N_t["Time (yr)"], N_t['Pu-241'], label='Pu-241')
ax.plot(N_t["Time (yr)"], N_t['Am-241'], label='Am-241')
ax.plot(N_fail["Time (yr)"], N_fail['moles'], ls='', marker='x', color ='k' ,label='canister fail')

ax.set_xlabel("time (years)") #, fontsize=9)
ax.set_ylabel("moles per cannister") #, fontsize=9)
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend()
fig.savefig('Output/Am_decay_chain_profiles.svg', transparent=False, bbox_inches="tight")

#%% plot of Am activity
yr_t_sec = 60 * 60 *24 *365 #seconds in a year
plt.style.use('Input/publish3.mplstyle')
fig, ax = plt.subplots(figsize=(5,5)) #size is in inches
ax.plot(N_t["Time (yr)"], N_t['Am-243'] * lamda[3] / yr_t_sec, label='Am-243')
# ax.plot(N_t["Time (yr)"], N_t['Cm-245'] * lamda[0] / yr_t_sec, label='Cm-245')
# ax.plot(N_t["Time (yr)"], N_t['Pu-241'] * lamda[1] / yr_t_sec, label='Pu-241')
ax.plot(N_t["Time (yr)"], N_t['Am-241'] * lamda[2] / yr_t_sec, label='Am-241')
ax.plot(N_t["Time (yr)"], N_t['Am-241'] * lamda[2] / yr_t_sec \
        + N_t['Am-243'] * lamda[3] / yr_t_sec, label='Total Am', color='grey')
ax.plot(N_fail.loc[2:3,"Time (yr)"], np.multiply(N_fail.loc[2:3,'moles'], [lamda[2],\
        lamda[3]]) / yr_t_sec, ls='', marker='x', color ='k' ,label='canister fail')

ax.set_xlabel("time (years)") #, fontsize=9)
ax.set_ylabel("activity per cannister (Becquerel)") #, fontsize=9)
# ax.set_yscale('log')
# ax.set_xscale('log')
ax.legend()
# fig.savefig('Output/Am_activity_profile.svg', transparent=False, bbox_inches="tight")

#%% function to round Round a number to a specified number of significant figures 
def round_to_sigfigs(x, sigfigs):
    if x == 0:
        return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (sigfigs - 1))

#%% calculate the concentration of the Am contaminated ground water

leach_fraction = assume.loc[0,'leach_fraction'] # fraction of  Am inventory leached per day upon canister failure
flow_rate = assume.loc[0, 'gw_flow'] # L/min

Am_lc_m = leach_fraction / 1440 / flow_rate * (N_fail.loc[2,'moles'] + N_fail.loc[3,'moles']) #Am concentration m/kg water
Darcy_flux = flow_rate * 0.001 / assume.loc[0,'rupture_area'] * 525960 # darcy flux through canister rupture in m/yr

#%% Make human readable report

text1 = f'The utilized ultrasonic evaluation of canister thickness was conducted {UT_df.loc[0,'t_past']} years after radionuclide inventory.\n \
HLW canister thickness remaining is {round((D*100),3)} cm.\n\
The current canister corrosion rate is {round(CorRate / 1e-6, 3)} μm/year.\n\
The canister is predicted to fail in {round((t_fail - UT_df.loc[0,'t_past']),1)} years ({round(t_fail, 1)} years after nuclide inventory).\n \n\
Tracked isotopes inventory at the time of predicted canister failure' 

text2 = f'Assuming Am fractional release rate of {leach_fraction} per day,\
and {flow_rate} L/min flow rate of groundwater through the failed canister,\
the concentration of  Am in the groundwater moving through the canister is {round_to_sigfigs(Am_lc_m,3)} m.'

report = 'Output/Canister_life_report.txt'
with open(report, mode='w', encoding='UTF-8') as output:
    print(text1, N_fail.loc[:,'Isotope':'moles'], ' ', text2, sep='\n', file=output)

#%% Update Crunchflow input file with Am concentration

file_source = 'Input/Am_RT_input_template.in'
file_destination = "RT_directory/RT_input.in"
pattern1 = "PythonValue1" #Key word for position of Am concentration in crunch flow template
pattern2 = "PythonValue2"
to_1 = str(round_to_sigfigs(Am_lc_m, 5))
to_2 = str(round_to_sigfigs(Darcy_flux, 5))
file = open(file_source,"r") # read template file
text = file.read()
file.close()

splitted_text = re.split(pattern1,text) #find patern and split sring at that point
modified_text1 = to_1.join(splitted_text)
splitted_text2 = re.split(pattern2,modified_text1)
modified_text = to_2.join(splitted_text2)
with open(file_destination, 'w') as file:
    file.write(modified_text)