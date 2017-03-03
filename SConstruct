"""
Build file for ebtel++
"""

import sys
import os

AddOption('--test',dest='test',action='store_true',help='Set this flag to compile the tests')
AddOption('--debug_compile',dest='debug_compile',action='store_true',help='Turn off optimizations for debugging.')
#allow user to specify custom path to libs and headers
AddOption('--libpath',dest='libpath',type='string',nargs=1,
        action='store',default=None,
        help='Comma-separated list of custom library paths if the defaults do not work.')
AddOption('--libs',dest='libs',type='string',nargs=1,
        action='store',default=None,
        help='Comma-separated list of custom libraries if the defaults do not work.')
AddOption('--includepath',dest='includepath',type='string',nargs=1,
        action='store',default=None,
        help='Comma-separated list of custom include paths if defaults do not work.')


subdirs = ['Radiation_Model','rsp_toolkit','source']
cxx_flags = ['-std=c++11']
try:
    CXX=os.environ['CXX']
except KeyError:
    CXX='g++'
if GetOption('debug_compile'):
    cxx_flags += ['-g','-Wall']
else:
    cxx_flags += ['-O3']
env = Environment(CXX=CXX,CXXFLAGS=cxx_flags)

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

if GetOption('libpath'):
    env['LIBPATH'] = GetOption('libpath').split(',')
if GetOption('includepath'):
    env['CPPPATH'] = GetOption('includepath').split(',')
if GetOption('libs'):
    env['LIBS'] = GetOption('libs').split(',')

allobjs = []
for sd in subdirs:
    consfile = os.path.join(sd,'SConscript')
    allobjs += env.SConscript(consfile,exports=['env'])

#TODO will this path always be resolved correctly?
if not os.path.exists('bin'):
    os.makedirs('bin')

print('Compiling ebtel++')
env.Program('bin/ebtel++.run',allobjs)

if GetOption('test'):
    print('Running tests...')
    try:
        python_path = os.path.join(os.environ['CONDA_PREFIX'],'bin','python')
    except KeyError:
        python_path = 'python'
    env.Command('dummy',None,'{} test/run_tests.py'.format(python_path))
