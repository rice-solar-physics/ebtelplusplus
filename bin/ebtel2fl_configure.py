#ebtel2fl_configure.py

#Will Barnes
#7 May 2015

#Import needed modules
import numpy as np
import os

class Configurer(object):
    
    def __init__(self,config_dictionary,root_dir,**kwargs):
        self.config_dictionary = config_dictionary
        self.root_dir = root_dir
        
        #Set up paths
        self.path_builder()
        
        if 'Hn' in kwargs:
            self.Hn = kwargs['Hn']
            
        if 'mc' in kwargs:
            self.mc = kwargs['mc']
        else:
            self.mc = False


    def path_builder(self,**kwargs):
        gen_path = self.config_dictionary['heat_species']+'_heating_runs/'
        if self.config_dictionary['amp_switch'] == 'uniform':
            gen_path = gen_path + 'alpha' + self.config_dictionary['amp_switch'] + '/'
        else:
            gen_path = gen_path + 'alpha' + str(np.fabs(self.config_dictionary['alpha'])) + '/'
        self.data_path = self.root_dir + gen_path + 'data/'
        self.config_path = self.root_dir + gen_path + 'config/'
        self.fn = 'ebtel2fl_L' + str(self.config_dictionary['loop_length']) + '_tn%d_tpulse' + str(2.0*self.config_dictionary['t_pulse_half']) + '_' + self.config_dictionary['solver']


    def print_xml_config(self,**kwargs):
        #Check if we have passed a filename
        #If not, pass a default filename
        if 'config_file' in kwargs:
            config_file = kwargs['config_file']
        else:
            config_file = '../config/ebtel_config.xml'

        #Open the file
        f = open(config_file,'w')

        #Print necessary header info
        f.write('<?xml version="1.0" ?>\n')
        f.write('<input>\n')

        #Loop through dictionary and print to xml file
        for key in self.config_dictionary:
            #Print tab delimiter, brackets and keyword
            f.write('\t<')
            f.write(key)
            f.write('>')
            #Check if entry is a list
            #If so print it as a list
            if isinstance(self.config_dictionary[key],list) or type(self.config_dictionary[key]).__name__ == 'ndarray':
                #Make temporary list
                temp = self.config_dictionary[key]
                #Skip to new line
                f.write('\n')
                #Begin loop over list
                for i in range(len(self.config_dictionary[key])):
                    f.write('\t\t<')
                    f.write(key+str(i))
                    f.write('>')
                    f.write(str(temp[i]))
                    f.write('</')
                    f.write(key+str(i))
                    f.write('>\n')
                #Print additional tab to preserve alignment
                f.write('\t')
            else:
                #Print value
                f.write(str(self.config_dictionary[key]))
            #Close the brackets and print newline
            f.write('</')
            f.write(key)
            f.write('>\n')

        #Close the main node of the file
        f.write('</input>')

        #Close the file
        f.close()


    def vary_wait_time(self,tn_a,tn_b,delta_tn,**kwargs):
        #Build wait time list
        t_wait = np.arange(tn_a,tn_b+delta_tn,delta_tn)
        #Iterate over wait times
        for i in range(len(t_wait)):
            self.stamp_arrays(t_wait[i])
            self.config_dictionary['h_nano'] = 2.0*self.Hn*self.config_dictionary['total_time']/(self.config_dictionary['num_events']*2.0*self.config_dictionary['t_pulse_half'])
            
            #Check if Monte-Carlo run
            if self.mc is not False:
                #Create directories in data and config if needed
                if not os.path.exists(self.config_path+self.fn%t_wait[i]):
                    os.makedirs(self.config_path+self.fn%t_wait[i])
                if not os.path.exists(self.data_path+self.fn%t_wait[i]):
                    os.makedirs(self.data_path+self.fn%t_wait[i])
                #Print config files for each mc run
                for j in range(self.calc_nmc):
                    self.config_dictionary['output_file'] = self.data_path+self.fn%t_wait[i]+'/'+self.fn%t_wait[i]+'_'+str(j)
                    self.print_xml_config(config_file=self.config_path+self.fn%t_wait[i]+'/'+self.fn%t_wait[i]+'_'+str(j)+'.xml')
            else:
                self.config_dictionary['output_file'] = self.data_path+self.fn%t_wait[i]
                self.print_xml_config(config_file=self.config_path+self.fn%t_wait[i]+'.xml')
                
                
    def calc_nmc(self,**kwargs):
        return int(np.ceil(self.mc/self.config_dictionary['num_events']))


    def stamp_arrays(self,ti,**kwargs):
        self.config_dictionary['num_events'] = int(np.ceil(self.config_dictionary['total_time']/(2.0*self.config_dictionary['t_pulse_half'] + ti)))
        self.config_dictionary['start_time_array'],self.config_dictionary['end_time_array'] = np.empty([self.config_dictionary['num_events']]),np.empty([self.config_dictionary['num_events']])
        if self.config_dictionary['amp_switch'] == 'file':
            self.config_dictionary['amp_array'] = np.empty([self.config_dictionary['num_events']])
        
        for i in range(self.config_dictionary['num_events']):
            self.config_dictionary['start_time_array'][i] = i*(2.0*self.config_dictionary['t_pulse_half'] + ti)
            self.config_dictionary['end_time_array'][i] = self.config_dictionary['start_time_array'][i] + 2.0*self.config_dictionary['t_pulse_half']
            if self.config_dictionary['amp_switch'] == 'file':
                x = np.random.rand(1)
                self.config_dictionary['amp_array'] = self.power_law_dist(x)


    def power_law_dist(self,x):
        return ((self.config_dictionary['amp1']**(self.config_dictionary['alpha']+1) - self.config_dictionary['amp0']**(self.config_dictionary['alpha']+1))*x + self.config_dictionary['amp0']**(self.config_dictionary['alpha']+1))**(1/(self.config_dictionary['alpha']+1))
        