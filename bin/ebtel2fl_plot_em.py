#ebtel2fl_plot_em.py

#Will Barnes
#14 May 2015

#Import needed modules
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit

class DEMPlotter(object):

    def __init__(self,temp_list,em_list,alpha,**kwargs):
        #static parameters
        self.fs = 18.0
        self.figsize = (10,10)
        self.linestyles = ('-','--','-.',':')
        self.format = 'eps'
        self.dpi = 1000
        #arguments
        self.temp_list = temp_list
        self.em_list = em_list
        self.alpha = alpha
        #keyword arguments
        if 'Tn' in kwargs:
            self.Tn = kwargs['Tn']
        else:
            self.Tn = np.arange(250,5250,250)
        self.Tndelta = self.Tn[1] - self.Tn[0]


    def plot_em_curves(self,**kwargs):
        #spacing between tn curves (artificial)
        delta_em = 0.2

        #set up figure
        fig = plt.figure(figsize=self.figsize)
        ax = fig.gca()

        #print lines
        for i in range(len(self.em_list)):
            if len(np.shape(np.array(self.em_list[i]))) > 1:
                mean_em = np.mean(self.inf_filter(self.em_list[i]),axis=0)
                temp_mean_em = np.array(mean_em)
                temp_mean_em[np.where(temp_mean_em==0.0)]=-np.float('Inf')
                mean_em = temp_mean_em
                mean_temp = np.mean(self.temp_list[i],axis=0)
                ax.plot(mean_temp,mean_em+i*delta_em,linestyle=self.linestyles[i%len(self.linestyles)],color='blue')
            else:
                ax.plot(np.array(self.temp_list[i]),np.array(self.em_list[i])+i*delta_em,linestyle=self.linestyles[i%len(self.linestyles)],color='blue')

        #set labels
        ax.set_title(r'EBTEL Two-fluid EM, $\alpha$ = '+str(self.alpha),fontsize=self.fs)
        ax.set_xlabel(r'$\log T$ (K)',fontsize=self.fs)
        ax.set_ylabel(r'$\log$EM (cm$^{-5}$)',fontsize=self.fs)
        ax.set_xlim([5.5,7.5])
        ax.set_ylim([27,33])

        #save or show the figure
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename']+'.'+self.format,format=self.format,dpi=self.dpi)
        else:
            plt.show()


    def plot_em_curve(self,tn_index,**kwargs):
        #get single list
        em_list = self.em_list[tn_index]
        temp_list = self.temp_list[tn_index]

        #set up figure
        fig = plt.figure(figsize=self.figsize)
        ax = fig.gca()

        #print lines
        mean_em = np.mean(em_list,axis=0)
        std_em = np.std(em_list,axis=0)
        mean_temp = np.mean(temp_list,axis=0)
        for i in range(len(temp_list)):
            ax.plot(temp_list[i],em_list[i],color='blue',linestyle=self.linestyles[-1])
        ax.plot(mean_temp,mean_em,color='black')
        ax.fill_between(mean_temp,mean_em-std_em,mean_em+std_em,facecolor='red',edgecolor='red',alpha=0.25)

        #set labels
        ax.set_title(r'EBTEL Two-fluid EM, $\alpha$ = '+str(self.alpha)+", $T_n$ = "+str(self.Tn[tn_index]),fontsize=self.fs)
        ax.set_xlabel(r'$\log T$ (K)',fontsize=self.fs)
        ax.set_ylabel(r'$\log$EM (cm$^{-5}$)',fontsize=self.fs)
        ax.set_xlim([5.5,7.5])
        ax.set_ylim([25,30])

        #save or show figure
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename']+'.'+self.format,format=self.format,dpi=self.dpi)
            plt.close('all')
        else:
            plt.show()


    def plot_em_max(self,temp_max,em_max,**kwargs):
        #set up figure
        fig = plt.figure(figsize=(self.figsize[0],0.7*self.figsize[1]))
        ax = fig.gca()
        ax_twin = ax.twinx()

        for i in range(len(self.Tn)):
            #calculate average and std
            mean_temp_max = np.mean(temp_max[i])
            std_temp_max = np.std(temp_max[i])
            mean_em_max = np.mean(em_max[i])
            std_em_max = np.std(em_max[i])
            #plot points
            ax.errorbar(self.Tn[i],mean_temp_max,yerr=std_temp_max,fmt='o',color='black')
            ax_twin.errorbar(self.Tn[i],mean_em_max,yerr=std_em_max,fmt='*',color='black')

        #set labels
        ax.set_title(r'EBTEL Two-fluid $T(\max(EM))$, $\alpha$ = '+str(self.alpha),fontsize=self.fs)
        ax.set_xlabel(r'$T_N$',fontsize=self.fs)
        ax.set_ylabel(r'$\log(T_{max})$',fontsize=self.fs)
        ax.set_ylim([5.5,7.0])
        ax.set_xlim([self.Tn[0]-self.Tndelta,self.Tn[-1]+self.Tndelta])
        ax_twin.set_ylabel(r'$\log$EM($T_{max}$) (cm$^{-5}$)',fontsize=self.fs)
        ax_twin.set_ylim([26,30])

        #save or show figure
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename']+'.'+self.format,format=self.format,dpi=self.dpi)
            plt.close('all')
        else:
            plt.show()


    def plot_em_slopes(self,a_cool,a_hot,**kwargs):
        #set up figure
        fig = plt.figure(figsize=self.figsize)
        ax = fig.gca()

        for i in range(len(self.Tn)):
            try:
                a_cool_mean = np.mean([a_cool[i][j] for j in np.where(np.array(a_cool[i]) != False)[0]])
                a_cool_std = np.std([a_cool[i][j] for j in np.where(np.array(a_cool[i]) != False)[0]])
                ax.errorbar(self.Tn[i],a_cool_mean,yerr=a_cool_std,fmt='o',color='blue')
            except:
                pass
            try:
                a_hot_mean = np.mean([a_hot[i][j] for j in np.where(np.array(a_hot[i]) != False)[0]])
                a_hot_std = np.std([a_hot[i][j] for j in np.where(np.array(a_hot[i]) != False)[0]])
                ax.errorbar(self.Tn[i],np.fabs(a_hot_mean),yerr=a_hot_std,fmt='o',color='red')
            except:
                pass

        #set labels
        ax.set_title(r'EBTEL Two-fluid Hot Shoulder Strength Comparison',fontsize=self.fs)
        ax.set_xlabel(r'$T_N$',fontsize=self.fs)
        ax.set_ylabel(r'$a_{hot,cool}$',fontsize=self.fs)
        ax.plot([self.Tn[0]-self.Tndelta,self.Tn[-1]+self.Tndelta],[2,2],'--k')
        ax.plot([self.Tn[0]-self.Tndelta,self.Tn[-1]+self.Tndelta],[3,3],'-k')
        ax.plot([self.Tn[0]-self.Tndelta,self.Tn[-1]+self.Tndelta],[5,5],'-.k')
        ax.set_ylim([0,10])
        ax.set_xlim([self.Tn[0]-self.Tndelta,self.Tn[-1]+self.Tndelta])

        #save or show figure
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename']+'.'+self.format,format=self.format,dpi=self.dpi)
            plt.close('all')
        else:
            plt.show()
            
            
    def inf_filter(self,nested_list,**kwargs):
        #preallocate space
        filtered_list = []
        #filter out infs in list and set to zero for averaging
        for i in nested_list:
            temp_array = np.array(i)
            temp_array[np.where(np.isinf(temp_array)==True)]=0.0
            filtered_list.append(temp_array)
        return filtered_list

