"""
Build file for ebtel++
"""

import sys
import os

subdirs = ['Radiation_Model','rsp_toolkit','source']

env = Environment(CXX='g++',CXXFLAGS=['-g','-Wall','-O3','-fno-stack-protector'])

if 'darwin' in sys.platform:
    print("Using Mac OS X compile options.")
    env.Append(CPPPATH=['/opt/local/include','/usr/include/malloc'])
    env.Append(LIBS=['boost_program_options-mt'])
    env.Append(LIBPATH=['/opt/local/lib'])
elif 'linux' in sys.platform:
    print("Using Linux compile options.")
    env.Append(CPPPATH=['/usr/include'])
    env.Append(LIBS=['boost_program_options'])
    env.Append(LIBPATH=['/usr/lib/x86_64-linux-gnu'])
else:
    print("Unrecognized platform. Using Windows compile options.")
    env.Append(CPPPATH=['/usr/local/include'])
    env.Append(LIBS=['boost_program_options'])
    env.Append(LIBPATH=['/usr/local/lib'])

allobjs = []
for sd in subdirs:
    consfile = os.path.join(sd,'SConscript')
    allobjs += env.SConscript(consfile,exports=['env'])

#TODO will this path always be resolved correctly?
if not os.path.exists('bin'):
    os.makedirs('bin')

env.Program('bin/ebtel++.run',allobjs)
