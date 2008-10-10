#!/bin/sh

#  Set up the build system -- need some symlinks to the build
#  directory.
#
if [ ! -e Makefile ] ; then
  if [ -e ../build/Makefile ] ; then
    ln -s ../build/Make.rules .
    ln -s ../build/Makefile .
  elif [ -e build/Makefile ] ; then
    ln -s build/Make.rules .
    ln -s build/Makefile .
  else
    echo "ERROR:  Couldn't find the Makefile!"
    exit 1
  fi
fi


#  If no target, try to figure out one based on uname.  This defaults to
#  the optimized target below.  If it works well, we can always use this
#  mechanism, and extend with "debug" or "profile" (e.g., "./configure.sh debug")
#
target=$1

if [ "x$target" = "xdebug" ] ; then
  opts="-debug";
  target=""
fi
if [ "x$target" = "xprofile" ] ; then
  opts="-profile";
  target=""
fi

if [ "x$target" = "x" ] ; then
  case `uname` in
    Darwin)
      target="Darwin-i386$opts"
      if [ "`uname -m`" = "Power Macintosh" ] ; then
          target="Darwin-ppc$opts"
      fi
      ;;
    FreeBSD)
      target="FreeBSD-i386$opts"
      if [ `uname -m` = "amd64" ] ; then
        target="FreeBSD-amd64$opts"
      fi
      ;;
    AIX)
      target="AIX$opts"
      ;;
    OSF1)
      target="OSF1$opts"
      ;;
    Linux)
      target="Linux-i686$opts"
      if [ `uname -m` = "x86_64" ] ; then
        target="Linux-amd64$opts"
      fi
      ;;
    *)
      echo "ERROR: Unknown uname of `uname` -- try manual configuration."
      exit 1
      ;;
  esac
fi


#
#  Look for the python headers.  We don't need the libraries.  This is
#  used by atac-driver/chainer only.
#

CFLAGS_PYTHON=

if [ x$PYTHONHOME != x ] ; then
  python=`echo $PYTHONHOME | cut -d: -f 1`
else
  python=`which python`
  if [ x$python != x ] ; then
    python=`dirname $python`
    python=`dirname $python`
  fi
fi

if [ x$python != x ] ; then
  if [ -e $python/include/python2.5/Python.h ]
  then
    CFLAGS_PYTHON="-I$python/include/python2.5"
  elif [ -e $python/include/python2.4/Python.h ]
  then
    CFLAGS_PYTHON="-I$python/include/python2.4"
  elif [ -e $python/include/python2.3/Python.h ]
  then
    CFLAGS_PYTHON="-I$python/include/python2.3"
  elif [ -e $python/include/python2.2/Python.h ]
  then
    CFLAGS_PYTHON="-I$python/include/python2.2"
  elif [ -e /usr/local/include/python2.5/Python.h ]
  then
    CFLAGS_PYTHON="-I/usr/local/include/python2.5"
  elif [ -e /usr/local/include/python2.4/Python.h ]
  then
    CFLAGS_PYTHON="-I/usr/local/include/python2.4"
  elif [ -e /usr/local/include/python2.3/Python.h ]
  then
    CFLAGS_PYTHON="-I/usr/local/include/python2.3"
  elif [ -e /usr/local/include/python2.2/Python.h ]
  then
    CFLAGS_PYTHON="-I/usr/local/include/python2.2"
  elif [ -e /System/Library/Frameworks/Python.framework/Versions/2.3/include/python2.3/Python.h ]
  then
    CFLAGS_PYTHON="-I/System/Library/Frameworks/Python.framework/Versions/2.3/include/python2.3"
  elif [ -e /usr/local/packages/python/include/python2.4/Python.h ]
  then
    CFLAGS_PYTHON="-I/usr/local/packages/python/include/python2.4"
  else
    echo "WARNING:  'python' program found, but can't find include file 'Python.h' -- is python installed correctly?"
    echo "WARNING:  Set PYTHONHOME to the directory containing include/Pyhton.h"
  fi
else
  echo "WARNING:  'python' program not found."
fi

if [ x$CFLAGS_PYTHON = x ] ; then
  echo "WARNING:  Will not build ATAC."
  WITHOUT_ATAC="atac-driver/ seatac/"
fi


case $target in
  Darwin-i386)
    echo "ERROR: Darwin Intel not supported because the developer doesn't have access to one."
    exit 1
    ;;
  Darwin-i386-debug)
    echo "ERROR: Darwin Intel not supported because the developer doesn't have access to one."
    exit 1
    ;;
  Darwin-ppc)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  OS-X, optimized
#
#  Using this usually generates internal compiler errors: -mpowerpc64
#    FIXED with the "November 2004" update!
#  Using this breaks dynamic library building:  -mdynamic-no-pic
#    We could have instead included -fPIC on the compile line.
#
# -malign-natural changes the size of structures compared to -g.
# -Wpadded supposedly warns when this happens.
#
FAST              := -fast -fPIC
FAST              := -O3 -funroll-loops -fstrict-aliasing -fsched-interblock -falign-loops=16 -falign-jumps=16 -falign-functions=16 -falign-jumps-max-skip=15 -falign-loops-max-skip=15 -malign-natural -ffast-math -mpowerpc-gpopt -force_cpusubtype_ALL -fstrict-aliasing -mtune=G5 -mcpu=G5
CC                := gcc
SHLIB_FLAGS       := -dynamiclib
CFLAGS_COMPILE    := \$(FAST) -fmessage-length=0 -D_THREAD_SAFE -Wall -Wno-char-subscripts
CLDFLAGS          := 
CLIBS             := 
CXX               := g++
CXXFLAGS_COMPILE  := \$(FAST) -fmessage-length=0 -D_THREAD_SAFE -Wall -Wno-char-subscripts
CXXLDFLAGS        := 
CXXLIBS           := 
CXXSHARED         := -Wl,-r -dynamic
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  Darwin-ppc-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  OS-X, debug
#
CC                := gcc
SHLIB_FLAGS       := -dynamiclib
CFLAGS_COMPILE    := -g3 -fmessage-length=0 -D_THREAD_SAFE -Wall -Wno-char-subscripts
CLDFLAGS          := 
CLIBS             := 
CXX               := g++
CXXFLAGS_COMPILE  := -g3 -fmessage-length=0 -D_THREAD_SAFE -Wall -Wno-char-subscripts
CXXLDFLAGS        := 
CXXLIBS           := 
CXXSHARED         := -dynamiclib
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  Darwin-ppc-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  OS-X, optimized, profiled (same as OS-X, optimized, but includes -pg)
#
FAST              := -pg -O3 -funroll-loops -fstrict-aliasing -fsched-interblock -falign-loops=16 -falign-jumps=16 -falign-functions=16 -falign-jumps-max-skip=15 -falign-loops-max-skip=15 -malign-natural -ffast-math -mpowerpc-gpopt -force_cpusubtype_ALL -fstrict-aliasing -mtune=G5 -mcpu=G5
CC                := gcc
SHLIB_FLAGS       := -dynamiclib
CFLAGS_COMPILE    := \$(FAST) -fmessage-length=0 -D_THREAD_SAFE -Wall -Wno-char-subscripts
CLDFLAGS          := -pg
CLIBS             := 
CXX               := g++
CXXFLAGS_COMPILE  := \$(FAST) -fmessage-length=0 -D_THREAD_SAFE -Wall -Wno-char-subscripts
CXXLDFLAGS        := -pg
CXXLIBS           := 
CXXSHARED         := -dynamiclib
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  FreeBSD-i386)
    rm -f Make.compilers
#    echo "Using linuxthreads by default!"
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized
THREADS           := -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads 
THREADL           := -llthread -llgcc_r
THREADS           := -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -O3 -fPIC \$(THREADS) -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL)
CXX               := g++
CXXFLAGS_COMPILE  := -O3 -fPIC \$(THREADS) -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL)
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  FreeBSD-i386-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized, profiled (same as FreeBSD, optimized, but includes -pg)
THREADS           := -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads 
THREADL           := -llthread -llgcc_r
THREADS           := -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -pg -O3 \$(THREADS) -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL)
CXX               := g++
CXXFLAGS_COMPILE  := -pg -O3 \$(THREADS) -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL)
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  FreeBSD-i386-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, debug, warnings
#  removed -Wredundant-decls cause it is annoying.
#  removed -Waggregate-return cause stl does it too much
#
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -g \$(THREADS) -fPIC -Wall -Wno-char-subscripts -Wshadow -Wpointer-arith -Wcast-qual \
  -Wcast-align -Wwrite-strings -Wconversion -Wstrict-prototypes -Wmissing-prototypes \
  -Wmissing-declarations -Wnested-externs  
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL)
CXX               := g++
CXXFLAGS_COMPILE  := -g \$(THREADS) -fPIC -Wall -Wno-char-subscripts -Wshadow -Wpointer-arith -Wcast-qual \
  -Wcast-align -Wwrite-strings -Wconversion
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL)
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  FreeBSD-amd64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized
THREADS           := -pthread
THREADL           := -pthread -lthr
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -O3 -fPIC \$(THREADS) -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL)
CXX               := g++
CXXFLAGS_COMPILE  := -O3 -fPIC \$(THREADS) -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL)
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  FreeBSD-amd64-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, debug, warnings
#  removed -Wredundant-decls cause it is annoying.
#  removed -Waggregate-return cause stl does it too much
#
THREADS           := -pthread
THREADL           := -pthread -lthr
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -g \$(THREADS) -fPIC -Wall -Wno-char-subscripts -Wshadow -Wpointer-arith -Wcast-qual \
  -Wcast-align -Wwrite-strings -Wconversion -Wstrict-prototypes -Wmissing-prototypes \
  -Wmissing-declarations -Wnested-externs  
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL)
CXX               := g++
CXXFLAGS_COMPILE  := -g \$(THREADS) -fPIC -Wall -Wno-char-subscripts -Wshadow -Wpointer-arith -Wcast-qual \
  -Wcast-align -Wwrite-strings -Wconversion
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL)
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  AIX)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, optimized
#  Tested, verified on 30jan03
CC                := xlc_r
SHLIB_FLAGS       := -G
CFLAGS_COMPILE    := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict -qcpluscmt
CLDFLAGS          := 
CLIBS             := 
CXX               := xlC_r
CXXFLAGS_COMPILE  := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict
CXXLDFLAGS        := 
CXXLIBS           := 
CXXSHARED         := -shared
ARFLAGS           := -X 64 ruv
INSTALL/          := $target/
EOF
    ;;
  AIX-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, optimized, profiled
#  Untested; might not work with optimization
CC                := xlc_r
SHLIB_FLAGS       := -G
CFLAGS_COMPILE    := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict -qcpluscmt
CLDFLAGS          := 
CLIBS             := 
CXX               := xlC_r
CXXFLAGS_COMPILE  := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict
CXXLDFLAGS        := 
CXXLIBS           := 
CXXSHARED         := -shared
ARFLAGS           := -X 64 ruv
INSTALL/          := $target/
EOF
    ;;
  AIX-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, debug
ccDBUG            := 
CC                := xlc_r
SHLIB_FLAGS       := -G
CFLAGS_COMPILE    := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath -qcpluscmt
CLDFLAGS          := 
CLIBS             := 
CXX               := xlC_r
CXXFLAGS_COMPILE  := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath
CXXLDFLAGS        := 
CXXLIBS           := 
CXXSHARED         := -shared
ARFLAGS           := -X 64 ruv
INSTALL/          := $target/
EOF
    ;;
  AIX-debug-full|aix-full-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, full debug, including memory checks
ccDBUG            := -qheapdebug -qinfo=uni:use -qcheck=all -qinitauto=FF -qflttrap=inv:ov:zero:en -qfloat=nans -DMEMORY_DEBUG
CC                := xlc_r
SHLIB_FLAGS       := -G
CFLAGS_COMPILE    := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath -qcpluscmt
CLDFLAGS          := \$(ccDBUG)
CLIBS             := 
CXX               := xlC_r
CXXFLAGS_COMPILE  := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath
CXXLDFLAGS        := \$(ccDBUG)
CXXLIBS           := 
CXXSHARED         := -shared
ARFLAGS           := -X 64 ruv
INSTALL/          := $target/
EOF
    ;;
  compaq|tru64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Tru64, native compilers, optimized
CC                := cxx
SHLIB_FLAGS       := -shared   # -G is something else
CFLAGS_COMPILE    := -D_REENTRANT -std -D_XOPEN_SOURCE=500 -D_OSF_SOURCE -pthread -w0 -fast
CDEPFLAGS         := -pthread -D_XOPEN_SOURCE=500
CLDFLAGS          := 
CLIBS             := -lpthread -lrt
CXX               := cxx
CXXFLAGS_COMPILE  := -D_REENTRANT -std -D_XOPEN_SOURCE=500 -D_OSF_SOURCE -pthread -w0 -fast
CXXDEPFLAGS       := -pthread -D_XOPEN_SOURCE=500
CXXLDFLAGS        := 
CXXLIBS           := -lpthread -lrt
CXXSHARED         := -shared
ARFLAGS           := ruv
CXX_TMP_ARCH      := cxx_repository
INSTALL/          := $target/
EOF
    ;;
  compaq-debug|tru64-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Tru64, native compilers, debug
CC                := cxx
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -D_REENTRANT -std -D_XOPEN_SOURCE=500 -D_OSF_SOURCE -pthread -g -w0 -trapuv
CLDFLAGS          := 
CLIBS             := -lpthread -lrt
CXX               := cxx
CXXFLAGS_COMPILE  := -D_REENTRANT -std -D_XOPEN_SOURCE=500 -D_OSF_SOURCE -pthread -gall -w0 -trapuv
CXXLDFLAGS        := 
CXXLIBS           := -lpthread -lrt
CXXSHARED         := -shared
ARFLAGS           := ruv
CXX_TMP_ARCH      := cxx_repository
INSTALL/          := $target/
EOF
    ;;
  Linux-i686)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux, optimized
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -O3 \$(THREADS) -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL) -ldl
CXX               := g++
CXXFLAGS_COMPILE  := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -O3 \$(THREADS) -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL) -ldl
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  Linux-amd64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux64, optimized
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -m64 -fPIC -D_REENTRANT -O3 \$(THREADS) -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL) -ldl
CXX               := g++
CXXFLAGS_COMPILE  := -m64 -fPIC -D_REENTRANT -O3 \$(THREADS) -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL) -ldl
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  Linux-amd64-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux64, optimized
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -m64 -fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g \$(THREADS) -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL) -ldl
CXX               := g++
CXXFLAGS_COMPILE  := -m64 -fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g \$(THREADS) -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL) -ldl
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  Linux-i686-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux, optimized
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g \$(THREADS) -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL) -ldl
CXX               := g++
CXXFLAGS_COMPILE  := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g \$(THREADS) -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL) -ldl
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;



  #  SUNLF needs to be set to allow for large file support on Solaris.  It
  #  should be whatever the following getconf's say.
  #
  #    getconf LFS_CFLAGS   = "-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64"
  #    getconf LFS_LDFLAGS  = ""
  #    getconf LFS_LFS_LIBS = ""
  #
  solaris)
    echo "Solaris is UNTESTED!"
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Solaris, gcc optimized
#
CC                := gcc
SHLIB_FLAGS       := -G #untested
CFLAGS_COMPILE    := -D_REENTRANT -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := 
CLIBS             := -lpthread -lrt
CXX               := g++
CXXFLAGS_COMPILE  := -D_REENTRANT -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer 
CXXLDFLAGS        := 
CXXLIBS           := -lpthread -lrt
CXXSHARED         := -shared
ARFLAGS           := ruv
INSTALL/          := $target/
EOF
    ;;
  *)
    echo "usage: $0 <configuration>"
    echo "          osx                       OS-X, optimized"
    echo "          osx-debug                 OS-X, debug"
    echo ""
    echo "          freebsd                   FreeBSD, optimized"
    echo "          freebsd-profile           FreeBSD, optimized, profiled"
    echo "          freebsd-debug             FreeBSD, debug and warnings"
    echo ""
    echo "          aix                       AIX, optimized"
    echo "          aix-profile               AIX, optimized, profiled (NOT TESTED)"
    echo "          aix-debug                 AIX, debug"
    echo ""
    echo "          tru64, compaq             Tru64, optimized"
    echo "          tru64-debug, compaq-debug Tru64, debug, warnings, trapuv"
    echo ""
    echo "          linux                     Linux, i686, optimized"
    echo "          linux64                   Linux, Opteron, optimized"
    echo ""
    echo "          solaris                   Solaris, gcc, optimized (STALE)"
    exit
    ;;
esac


cat <<EOF >> Make.compilers
PERL              := /usr/bin/perl
.EXE              := 
.SO               := .so
.A                := .a
.O                := .o
CLD               := \${CC}
CXXLD             := \${CXX}
CCDEP		  := gcc -MM -MG
CXXDEP	          := g++ -MM -MG
CLIBS             += -lm -lbz2
CXXLIBS           += -lm -lbz2
CFLAGS_PYTHON     := $CFLAGS_PYTHON
WITHOUT           := $WITHOUT_ATAC
EOF

if [ ! -e ./makepath ] ; then
  echo "include Make.compilers"       >> makepath-makefile
  echo "makepath: build/makepath.c"   >> makepath-makefile
  echo "	\${CC} \${CFLAGS} \${CLDFLAGS} -o ./makepath build/makepath.c \${CLIBS}" >> makepath-makefile
  gmake -f makepath-makefile
  rm -f makepath-makefile
fi

echo "Configured."

#cat Make.compilers
