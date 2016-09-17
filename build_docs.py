"""
ebtel++ docs build script using cldocs
"""

import glob
import os
import subprocess

#Get all headers and sources
sources=glob.glob('source/*.cpp')
headers=glob.glob('source/*.h')

#check for travis
on_travis = os.environ.get('TRAVIS',None)

#Set build options
if on_travis is not None and on_travis.lower()=='true':
    build_opts = '--report --merge docs/ext --output docs/html '
else:
    build_opts = '--static --report --merge docs/ext --output docs/html '
build_opts += ' '.join(sources) + ' ' + ' '.join(headers)

#Set needed CXX flags
cxx_flags='-I/opt/local/include -I/usr/include/malloc'

#set path to cldoc exec (if not already in path)
exec_path='cldoc'

#Run command
subprocess.call(exec_path+' generate ' + cxx_flags + ' -- ' + build_opts,shell=True)
