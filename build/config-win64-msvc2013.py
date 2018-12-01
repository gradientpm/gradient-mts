

BUILDDIR        = '#build/release'
DISTDIR         = '#dist'
CXX             = 'cl'
CC              = 'cl'
# /O2=optimize for speed, global optimizations, intrinsic functions, favor fast code, frame pointer omission
# /EHsc=C++ exceptions, /fp:fast=Enable reasonable FP optimizations, /GS-=No buffer security checks, /GL=whole program optimizations
# To include debug information add '/Z7' to CXXFLAGS and '/DEBUG' to LINKFLAGS
CXXFLAGS        = ['/nologo', '/O2', '/fp:fast', '/D', 'WIN32', '/D', 'WIN64', '/W3', '/EHsc', '/GS-', '/GL', '/MD', '/D', 'MTS_DEBUG', '/D', 'DOUBLE_PRECISION', '/D', 'SPECTRUM_SAMPLES=3', '/D',  '_CONSOLE', '/D', 'NDEBUG', '/D', 'OPENEXR_DLL', '/openmp', '/D', 'MTS_OPENMP']
SHCXXFLAGS      = CXXFLAGS
TARGET_ARCH     = 'x86_64'
MSVC_VERSION    = '12.0'
LINKFLAGS       = ['/nologo', '/SUBSYSTEM:CONSOLE', '/MACHINE:X64', '/FIXED:NO', '/OPT:REF', '/OPT:ICF', '/LTCG', '/NODEFAULTLIB:LIBCMT', '/MANIFEST']
BASEINCLUDE     = ['#include', '#dependencies/vs2013/include']
BASELIB         = ['msvcrt', 'ws2_32', 'Half']
BASELIBDIR      = ['#dependencies/vs2013/lib/x64_vc12']
OEXRINCLUDE     = ['#dependencies/vs2013/include/openexr']
OEXRLIB         = ['IlmImf', 'IlmThread', 'Iex', 'zlib']
BOOSTLIB        = ['boost_system-vc120-mt-1_53', 'boost_filesystem-vc120-mt-1_53', 'boost_thread-vc120-mt-1_53']
COLLADAINCLUDE  = ['#dependencies/vs2013/include/collada-dom', '#dependencies/vs2013/include/collada-dom/1.4']
COLLADALIB      = ['libcollada14dom24']
XERCESLIB       = ['xerces-c_3']
PNGLIB          = ['libpng16']
JPEGLIB         = ['jpeg']
GLLIB           = ['opengl32', 'glu32', 'glew32mx', 'gdi32', 'user32']
GLFLAGS         = ['/D', 'GLEW_MX']
SHLIBPREFIX     = 'lib'
SHLIBSUFFIX     = '.dll'
LIBSUFFIX       = '.lib'
PROGSUFFIX      = '.exe'
PYTHON27LIB     = ['boost_python27-vc120-mt-1_53', 'python27']
PYTHON27INCLUDE = ['#dependencies/vs2013/include/python27']
PYTHON32LIB     = ['boost_python32-vc120-mt-1_53', 'python32']
PYTHON32INCLUDE = ['#dependencies/vs2013/include/python32']
PYTHON33LIB     = ['boost_python33-vc120-mt-1_53', 'python33']
PYTHON33INCLUDE = ['#dependencies/vs2013/include/python33']
PYTHON34LIB     = ['boost_python34-vc120-mt-1_53', 'python34']
PYTHON34INCLUDE = ['#dependencies/vs2013/include/python34']
QTINCLUDE       = ['#dependencies/vs2013/qt/include']
QTDIR           = '#dependencies/vs2013/qt/x64_vc12'
FFTWLIB         = ['libfftw-3.3']
