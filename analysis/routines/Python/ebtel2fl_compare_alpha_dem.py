import numpy as np
import matplotlib.pyplot as plt

#Set parameters
L = 120.0
tpulse = 200.0
solver = 'rka4'
species = 'electron'

#Set the root directory
root_dir = '/data/datadrive2/EBTEL-2fluid_runs/'+species+ '_heating_runs/'

#Declare alpha vector
alpha = [1.5,2.0,2.5]
colors = ['r','b','g']

#Set up the figure
fig = plt.figure(figsize=(10,10))
ax = fig.gca()
fs = 18

#Loop over alpha values
for i in range(len(alpha)):
    #Load the data
    matlong = np.loadtxt(root_dir+'alpha'+str(alpha[i])+'/'+'data/ebtel2fl_L'+str(L)+'_tn5000_tpulse200.0_'+solver+'_dem.txt')
    matshort = np.loadtxt(root_dir+'alpha'+str(alpha[i])+'/'+'data/ebtel2fl_L'+str(L)+'_tn250_tpulse200.0_'+solver+'_dem.txt')
    #Do the plotting
    ax.plot(matlong[:,0],matlong[:,2],color=colors[i],linestyle='-')
    ax.plot(matshort[:,0],matshort[:,2],color=colors[i],linestyle='-.')
    
ax.set_title(r'L = '+str(L)+' Mm',fontsize=fs)
ax.set_xlabel(r'$\log(T_{DEM})$ (K)',fontsize=fs)
ax.set_ylabel(r'$\log($DEM$)$ (cm$^{-5}$ K$^{-1}$)',fontsize=fs)
ax.set_xlim([5.5,7.5])

#Show the plot
plt.show()
