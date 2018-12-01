import os, sys

BUILDDIR       = '#build/release'
DISTDIR        = '#dist'
CXX            = 'g++'
CC             = 'gcc'
CXXFLAGS       = ['-O3', '-Wall', '-g', '-pipe', '-march=nocona', '-msse4.2', '-ftree-vectorize', '-mfpmath=sse', '-funsafe-math-optimizations', '-fno-rounding-math', '-fno-signaling-nans', '-fno-math-errno', '-fomit-frame-pointer', '-DMTS_DEBUG', '-DSINGLE_PRECISION', '-DSPECTRUM_SAMPLES=3', '-DMTS_SSE', '-DMTS_HAS_COHERENT_RT', '-fopenmp', '-fvisibility=hidden', '-mtls-dialect=gnu2', '-std=gnu++11']
LINKFLAGS      = []
SHLINKFLAGS    = ['-rdynamic', '-shared', '-fPIC', '-lstdc++']
BASEINCLUDE    = ['#include']
BASELIB        = ['dl', 'm', 'pthread', 'gomp']
BASELIBDIR     = ['/home/neodym60/renderers/lib']
EIGENINCLUDE   = ['/software/CentOS-6/eb/software/Core/Eigen/3.2.8/include']
OEXRINCLUDE    = ['/home/neodym60/renderer/include/OpenEXR']
OEXRLIB        = ['Half', 'IlmImf', 'z']
PNGLIB         = ['png']
JPEGLIB        = ['jpeg']
XERCESINCLUDE  = ['/software/CentOS-6/eb/software/Toolchain/iomkl/2015b/Xerces-C++/3.1.4/include']
XERCESLIBDIR   = '/software/CentOS-6/eb/software/Toolchain/iomkl/2015b/Xerces-C++/3.1.4/lib'
XERCESLIB      = ['xerces-c']
GLLIB          = ['GL', 'GLU', 'GLEWmx', 'Xxf86vm', 'X11']
GLFLAGS        = ['-DGLEW_MX']
BOOSTINCLUDE   = ['/home/neodym60/renderers/include']
BOOSTLIBDIR    = ['/home/neodym60/renderers/lib']
BOOSTLIB       = ['boost_system', 'boost_filesystem', 'boost_thread']
GLINCLUDE      = '/home/neodym60/renderers/include'
#COLLADAINCLUDE = ['/usr/include/collada-dom', '/usr/include/collada-dom/1.4']
#COLLADALIB     = ['collada14dom', 'xml2']
FFTWLIB        = ['fftw3_threads', 'fftw3']

PYTHON26LIB     = ['boost_python', 'python2.6']
PYTHON26INCLUDE = ['/usr/include/python2.6']

# The following runs a helper script to search for installed Python
# packages that have a Boost Python library of matching version.
# A Mitsuba binding library will be compiled for each such pair.
# Alternatively, you could also specify the paths and libraries manually
# using the variables PYTHON27INCLUDE, PYTHON27LIB, PYTHON27LIBDIR etc.

import sys, os
sys.path.append(os.path.abspath('../data/scons'))
from detect_python import detect_python
locals().update(detect_python())

