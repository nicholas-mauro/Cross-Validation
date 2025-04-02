# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:41:15 2023

@author: maurna
"""

# Import Data from File into array and plots the data.  Data from Matlab.  No headers.  Time (ps), delta N
Rootfile = '20mev_3.75.txt'
infile = open(Rootfile) #Open the file 
outTitle = Rootfile +'_no_line.txt' #Strip off the blank lines
target = open(outTitle, 'w')
print(outTitle)
for line in infile:
    if not line.strip():
        continue
    else:
        print(line)
        target.write(line)

print("And finally, we close it.")
target.close()

Time_list = []  # list of TIMES
DN_list = []  #List of Delta N
adata=[]
taudata=[]  #setting up arrays for analysis
betadata=[]
tau1data=[]
tau2data=[]
a1data=[]
a2data=[]
line_count = 0
dummy = 0
target = open(outTitle, 'r')

for line in target:  #reading the time (dummy[0]) and the RK metric (dummy [1])
   dummy = line.split()
   print(dummy)
   Time_list.append(float(dummy[0]))
   DN_list.append(float(dummy[1]))
target.close()


#check the lengths of the list
len(Time_list)==len(DN_list)

import matplotlib.pyplot as plt
plt.plot(Time_list,DN_list)


import random
import math

#input length of array, fraction in the training set

length = len(Time_list)

def Training(length,fraction):  #this function takes as imput the length of the data and fraction that is wanted and outputs a random list of indices for training data
    a=[]
    afirst=[]
    alast=[]
    for i in range(0,length):
        a.append(i)
    random.shuffle(a)
    for i in range(0,math.ceil(fraction*len(a))):
        afirst.append(a[i])
    for i in range(math.ceil(fraction*len(a)),len(a)):
        alast.append(a[i])
    afirst.sort()
    alast.sort()
    return afirst, alast

trial1,trial2=Training(length, .8)  #outputs example list of training,cross-check

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def stretched(t, a, tau, beta):     #stretched exponential with no offset
    return a * np.exp((-1)*(t/tau)**beta)

def twoexp(t,a1,tau1,a2,tau2):
    return a1*np.exp(-t/tau1)+a2*np.exp(-t/tau2) #two exponentials with no offset

def SelectData(Time_list,DN_list,trial1,trial2):
    Training_time=[]
    Training_DN = []
    Cross_time=[]
    Cross_DN = []
    for i in trial1:
        Training_time.append(Time_list[i])
        Training_DN.append(DN_list[i])
    for i in trial2:
        Cross_time.append(Time_list[i])
        Cross_DN.append(DN_list[i])
    return Training_time,Training_DN,Cross_time,Cross_DN  #output with be four arrays, extraning training/cross-check

Data = SelectData(Time_list, DN_list, trial1, trial2)  #output with be four arrays, extraning training/cross-check

"""
Outputs here are SSR.  Need different function for each fitting

"""


def Sum_Squares_Residuals_stretched(time_array,Experimental_Data,function,popt):  #using minimizing the sum of the squared of the residuals for finding best fit.
    sum1=0
    for i in range(0,len(time_array)):
        a=(Experimental_Data[i]-function(time_array[i],popt[0],popt[1],popt[2]))**2
        sum1=sum1+(Experimental_Data[i]-function(time_array[i],popt[0],popt[1],popt[2]))**2
        #print(a)
    return sum1
def Sum_Squares_Residuals_twoexp(time_array,Experimental_Data,function,popt1):
    sum1=0
    for i in range(0,len(time_array)):
        a=(Experimental_Data[i]-function(time_array[i],popt1[0],popt1[1],popt1[2],popt1[3]))**2
        sum1=sum1+(Experimental_Data[i]-function(time_array[i],popt1[0],popt1[1],popt1[2],popt1[3]))**2
        #print(a)
    return sum1

Number_trials = 500  #chose some large number of trials.  Each trial randomly picks training data
fraction = 0.8

matrix = [[None for x in range(13)]  # 11 columns: Trial, Training Time, Training DN, Cross_Time, Cross DN, Popt, Popt2, SSR1, SSR2
             for y in range(Number_trials+1)] 
matrix[0] = ['Trial', 'Training Time', 'Training DN', 'Cross_Time', 'Cross DN', 'Popt', 'Popt1', 'SSR_training', 'SSR1_training','SSR_cross', 'SSR1_cross', 'SSR (Stretched) Full', 'SSR1 (Two Exponential) Full'] #header
"""

Matrix elements may be arrays themselves.  Data can be output to a file as necessary
popt is the stretched
popt1 is the two

"""
for i in range (1,Number_trials+1):
    matrix[i][0]=i
    trial1,trial2=Training(length, .8) #get shuffled or random breakdown of training/cross data
    Data = SelectData(Time_list, DN_list, trial1, trial2) #output time training, DN training, Time cross, DN cross
    matrix[i][1]=Data[0]
    matrix[i][2]=Data[1]
    matrix[i][3]=Data[2]
    matrix[i][4]=Data[3]
    popt, pcov = curve_fit(stretched, Data[0], Data[1],p0=[3,1,1])
    popt1, pcov1 = curve_fit(twoexp, Data[0], Data[1],p0=[2.5, 1.3, 0.5,1.3], bounds=(0,[100,100,100,100]))
    taudata.append(popt[1])
    adata.append(popt[0])
    betadata.append(popt[2])
    if popt1[1]<=popt1[3]:
        tau1data.append(popt1[1])
        tau2data.append(popt1[3])
        a1data.append(popt1[0])
        a2data.append(popt1[2])
    else:
        tau1data.append(popt1[3])
        tau2data.append(popt1[1])
        a1data.append(popt1[2])
        a2data.append(popt1[0])
    matrix[i][5]=popt
    matrix[i][6]=popt1
    matrix[i][7] = Sum_Squares_Residuals_stretched(Data[0], Data[1], stretched, popt)
    matrix[i][8] = Sum_Squares_Residuals_twoexp(Data[0], Data[1], twoexp, popt1)
    matrix[i][9] = Sum_Squares_Residuals_stretched(Data[2], Data[3], stretched, popt)
    matrix[i][10] = Sum_Squares_Residuals_twoexp(Data[2], Data[3], twoexp, popt1)
    matrix[i][11] = Sum_Squares_Residuals_stretched(Time_list,DN_list, stretched,popt)
    matrix[i][12] = Sum_Squares_Residuals_twoexp(Time_list, DN_list, twoexp, popt1)

"""
for i in range(1,Number_trials+1):
    print(matrix[i][5]) 

for i in range(1,Number_trials+1):
    print(matrix[i][6]) 
"""

SSR_training = 0 
for i in range(1,Number_trials+1):
    print(matrix[i][7])
    
    SSR_training = SSR_training+matrix[i][7]     
print('ok')
SSR1_training = 0
for i in range(1,Number_trials+1):
    print(matrix[i][8])
    
    SSR1_training = SSR1_training+matrix[i][8]
print('ok')
SSR_cross = 0
for i in range(1,Number_trials+1):
    print(matrix[i][9])
    SSR_cross = SSR_cross+matrix[i][9]     
print('ok')
SSR1_cross = 0
for i in range(1,Number_trials+1):
    print(matrix[i][10])
    SSR1_cross = SSR1_cross+matrix[i][10]

SSR_Total = 0
for i in range(1,Number_trials+1):
    print(matrix[i][11])
    SSR_Total = SSR_Total+matrix[i][11]

SSR1_Total = 0
for i in range(1,Number_trials+1):
    print(matrix[i][12])
    SSR1_Total = SSR1_Total+matrix[i][12]

print('Stretched Training,  Two Exponential Training,  Stretched Cross Check,   Two Exponential Cross Check,  Stretched total,  Two Exponential Total, Ave Tau, STDEV Tau, Ave Tau1, STDEV Tau1, Ave Tau2, STDEV Tau2, a, STDEV a, a1, STDEV a1, a2, STDEV a2')
print(SSR_training,SSR1_training,SSR_cross,SSR1_cross, SSR_Total, SSR1_Total) #Sum of th squares (SSR), SSR = streched, SSR1 = two exp, training is the subset of data, ssr1 is the rest of the data
Test = SSR_Total-SSR_training-SSR_cross
Test1 = SSR1_Total-SSR1_training-SSR1_cross

#print('Differences in SSR and SSR1')
#print(Test, Test1)
print(popt,popt1)
print(SSR_training,SSR1_training,SSR_cross,SSR1_cross, SSR_Total, SSR1_Total,np.average(taudata),np.std(taudata),np.average(tau1data),np.std(tau1data),np.average(tau2data),np.std(tau2data),np.average(adata),np.std(adata),np.average(a1data),np.std(a1data),np.average(a2data),np.std(a2data))
plt.plot(Data[0],Data[1],marker="+",label = 'Fitting Data', markersize = 10)
#plt.plot(Data[2],Data[3],'b+', label = 'Training Data')
plt.plot(Data[2],Data[3],marker=".", label = 'Training Data',markersize = 20)         
#plt.plot(Time_list,stretched(Time_list, popt[0],popt[1],popt[2]), label = "Stretched Exponential Fit")
Twoexp_data = []
for i in range(0,length):
    #print(i)
    Twoexp_data.append(twoexp(Time_list[i], popt1[0], popt1[1], popt1[2], popt1[3]))
    
#plt.plot(Time_list,Twoexp_data, label = "Two Exponential Fit")

"""
Sum_Squares_Residuals_stretched(Data[0], Data[1], stretched, popt)
Sum_Squares_Residuals_twoexp(Data[0], Data[1], twoexp, popt1)
"""
#checking SSR calulcation manually
#popt, pcov = curve_fit(stretched, xdata, ydata,p0=[3,1,1,1])
#popt1, pcov1 = curve_fit(twoexp, xdata, ydata,p0=[2.5, 1.3, 0.5,1.3, .5])

plt.title('RK Metric Example, Starting Fit at 3.75 1/Angstroms')
plt.ylabel('RK Metric')
plt.xlabel('Time (Picoseconds)')
leg = plt.legend()
#construct training set and target set