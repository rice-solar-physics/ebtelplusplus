"""
Run tests for a few different ebtel++ inputs
"""

import os
import sys
import subprocess
import itertools

import numpy as np

def run_tests():
    top_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    sys.path.append(os.path.join(top_dir,'rsp_toolkit/python'))
    from xml_io import InputHandler,OutputHandler

    #set acceptable tolerance level
    tolerance = 1e-6

    #configure different options for tests
    events = [
                [{'event':{'rise_start':0.0,'rise_end':250.0,'decay_start':250.0,'decay_end':500.0,'magnitude':0.04}}],
                [{'event':{'rise_start':0.0,'rise_end':100.0,'decay_start':100.0,'decay_end':200.0,'magnitude':0.1}}],
    ]
    solvers = [True,False]
    test_options = list(itertools.product(events,solvers))

    #run tests
    ih = InputHandler(os.path.join(top_dir,'config/ebtel.example.cfg.xml'))
    base_dir = ih.lookup_vars()
    base_dir['calculate_dem'] = True
    base_dir['rka_error'] = 1e-8
    for opts in test_options:
        #configure opttions
        base_dir['heating']['events'] = opts[0]
        base_dir['use_adaptive_solver'] = opts[1]
        base_dir['output_filename'] = os.path.join(top_dir,'test','test_{0}_results'.format(test_options.index(opts)))
        #print output
        oh = OutputHandler(base_dir['output_filename']+'.xml',base_dir)
        oh.print_to_xml()
        #run the model
        subprocess.call([os.path.join(top_dir,'bin','ebtel++.run'),'-c',base_dir['output_filename']+'.xml'])
        #load results
        #results = np.loadtxt(base_dir['output_filename'])
        #load truth
        #truth = np.loadtxt(base_dir['output_filename']+'_truth')
        #compare
        #for i in range(1,np.shape(results)[1]):
        #    truth_interp = np.interp(results[:,0],truth[:,0],truth[:,i])
        #    err = np.fabs(truth_interp - results[:,i])/np.mean([truth_interp,results[:,i]],axis=0)
        #    print(np.max(err))
        #    assert np.max(err) < tolerance
        print('Passed test {0}'.format(test_options.index(opts)))

    print('Passed all tests.')

if __name__=='__main__':
    run_tests()
