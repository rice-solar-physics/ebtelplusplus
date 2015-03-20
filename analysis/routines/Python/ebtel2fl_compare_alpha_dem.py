import numpy as np
import matplotlib.pyplot as plt

#Set parameters
L = 120.0
tpulse = 100.0
solver = 'rka4'
species = 'ion'
tnlong = 5000
tnshort = 250

#Set the root directory
root_dir = '/data/datadrive2/EBTEL-2fluid_runs/'+species+ '_heating_runs/'

#Declare alpha vector
alpha = [1.5,2.0,2.5]
colors = ['r','b','g']
lines = []

#Set up the figure
fig = plt.figure(figsize=(10,10))
ax = fig.gca()
fs = 18

#Loop over alpha values
for i in range(len(alpha)):
    #Load the data
    matlong = np.loadtxt(root_dir+'alpha'+str(alpha[i])+'/'+'data/ebtel2fl_L'+str(L)+'_tn' + str(tnlong) + '_tpulse'+str(tpulse)+'_'+solver+'_dem.txt')
    matshort = np.loadtxt(root_dir+'alpha'+str(alpha[i])+'/'+'data/ebtel2fl_L'+str(L)+'_tn' + str(tnshort) + '_tpulse'+str(tpulse)+'_'+solver+'_dem.txt')
    #Do the plotting
    #Make sure the color correspondence is correct
    line = ax.plot(matlong[:,0],matlong[:,4],color=colors[i],linestyle='-',label=r'$\alpha$ = '+str(alpha[i]))
    ax.plot(matshort[:,0],matshort[:,4],color=colors[i],linestyle='-.')
    lines = line + lines

ax.set_title(r'L = '+str(L)+' Mm',fontsize=fs)
ax.set_xlabel(r'$\log(T)$ (K)',fontsize=fs)
ax.set_ylabel(r'$\log($EM$)$ (cm$^{-5}$)',fontsize=fs)
ax.set_xlim([5.5,7.5])
labels = [l.get_label() for l in lines]
ax.legend(lines,labels,loc=4)

#Show the plot
plt.show()
