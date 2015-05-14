#ebtel2fl_plot_em.py

#Will Barnes
#14 May 2015

#Import needed modules
import numpy as np
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
        self.dpi = '1000'
        #arguments
        self.temp_list = temp_list
        self.em_list = em_list
        self.alpha = alpha
        #keyword arguments
        if 'Tn' not in kwargs:
            self.Tn = np.arange(250,5250,250)
        else:
            self.Tn = kwargs['Tn']
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
                mean_em = np.mean(self.em_list[i],axis=0)
                std_em = np.std(self.em_list[i],axis=0)
                mean_temp = np.mean(self.temp_list[i],axis=0)
                ax.plot(mean_temp,mean_em+i*delta_em,linestyle=self.linestyles[i%len(self.linestyles)],color='blue')
                #ax.fill_between(mean_temp,mean_em+i*delta_em-std_em,mean_em+i*delta_em+std_em,facecolor='red',alpha=0.25)
            else:
                ax.plot(self.temp_list[i],self.em_list[i]+i*delta_em,linestyle=self.linestyles[i%len(self.linestyles)],color='blue')

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
        ax.fill_between(mean_temp,mean_em-std_em,mean_em+std_em,facecolor='red',alpha=0.25)

        #set labels
        ax.set_title(r'EBTEL Two-fluid EM, $\alpha$ = '+str(self.alpha),fontsize=self.fs)
        ax.set_xlabel(r'$\log T$ (K)',fontsize=self.fs)
        ax.set_ylabel(r'$\log$EM (cm$^{-5}$)',fontsize=self.fs)
        ax.set_xlim([5.5,7.5])
        ax.set_ylim([24,29])

        #save or show figure
        if 'print_fig_filename' in kwargs:
            plt.savefig(kwargs['print_fig_filename']+'.'+self.format,format=self.format,dpi=self.dpi)
        else:
            plt.show()


    def plot_em_max(self,**kwargs):
        #set up figure
        fig = plt.figure(figsize=(self.figsize[0],0.7*self.figsize[1]))
        ax = fig.gca()
        ax_twin = ax.twinx()

        for i in range(len(self.Tn)):
            temp_max = []
            em_max = []
            #calculate max values
            for j in range(len(self.em_list[i])):
                i_max = np.argmax(self.em_list[i][j])
                temp_max.append(self.temp_list[i][j][i_max])
                em_max.append(self.em_list[i][j][i_max])
            #calculate average and std
            mean_temp_max = np.mean(temp_max)
            std_temp_max = np.std(temp_max)
            mean_em_max = np.mean(em_max)
            std_em_max = np.std(em_max)
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
        else:
            plt.show()


    def plot_em_slopes(self,a_cool,a_hot,**kwargs):
        #set up figure
        fig = plt.figure(figsize=(1.2*self.figsize[0],self.figsize[1]))
        ax = fig.gca()

        for i in range(len(self.Tn)):
            #try:
                a_cool_mean = np.mean([a_cool[i][j] for j in np.where(np.array(a_cool[i])[0] != False)])
                a_cool_std = np.std([a_cool[i][j] for j in np.where(np.array(a_cool[i])[0] != False)])
                ax.errorbar(self.Tn[i],a_cool_mean,yerr=a_cool_std,fmt='o',color='blue')
            #except:
                #pass
            #try:
                a_hot_mean = np.mean([a_hot[i][j] for j in np.where(np.array(a_hot[i])[0] != False)])
                a_hot_std = np.std([a_hot[i][j] for j in np.where(np.array(a_hot[i])[0] != False)])
                ax.errorbar(self.Tn[i],a_hot_mean,yerr=a_hot_std,fmt='o',color='red')
            #except:
                #pass

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
        else:
            plt.show()
