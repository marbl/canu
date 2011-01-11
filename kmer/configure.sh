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
      if [ `uname -m` = "x86_64" ] ; then
        target="Darwin-amd64$opts"
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
      if [ `uname -m` = "ia64" ] ; then
        target="Linux-ia64$opts"
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
  if [ -e $python/include/python2.6/Python.h ]
  then
    CFLAGS_PYTHON="-I$python/include/python2.6"
  elif [ -e $python/include/python2.5/Python.h ]
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
  elif [ -e /usr/local/include/python2.6/Python.h ]
  then
    CFLAGS_PYTHON="-I/usr/local/include/python2.6"
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

#
#  Decide on compilers to use.  Unfortunately, all the options are tuned for gcc/g++.
#  In particular, -m64 and -W* and -f* aren't liked by Intel compilers.
#

if [ x$CC = x ] ; then
    CC="gcc"
fi

if [ x$CXX = x ] ; then
    CXX="g++"
fi

#
#  Emit architecture specific configurations.
#

case $target in
  Darwin-i386|Darwin-amd64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  OS-X, optimized
#
CC                := $CC
SHLIB_FLAGS       := -dynamiclib
CFLAGS_COMPILE    := -fast -fPIC -m64 -fmessage-length=0 -D_REENTRANT -D_THREAD_SAFE -Wall -Wno-char-subscripts
CLDFLAGS          := -m64
CLIBS             := 
CXX               := $CXX
CXXFLAGS_COMPILE  := -fast -fPIC -m64 -fmessage-length=0 -D_REENTRANT -D_THREAD_SAFE -Wall -Wno-char-subscripts
CXXLDFLAGS        := -m64
CXXLIBS           := 
CXXSHARED         := -Wl,-r -dynamic
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  Darwin-i386-debug|Darwin-amd64-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  OS-X, debug
#
CC                := $CC
SHLIB_FLAGS       := -dynamiclib
CFLAGS_COMPILE    := -g3 -m64 -fmessage-length=0 -D_REENTRANT -D_THREAD_SAFE -Wall -Wno-char-subscripts
CLDFLAGS          := -m64
CLIBS             := 
CXX               := $CXX
CXXFLAGS_COMPILE  := -g3 -m64 -fmessage-length=0 -D_REENTRANT -D_THREAD_SAFE -Wall -Wno-char-subscripts
CXXLDFLAGS        := -m64
CXXLIBS           := 
CXXSHARED         := -Wl,-r -dynamic
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
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
CC                := $CC
SHLIB_FLAGS       := -dynamiclib
CFLAGS_COMPILE    := \$(FAST) -fmessage-length=0 -D_REENTRANT -D_THREAD_SAFE -Wall -Wno-char-subscripts
CLDFLAGS          := 
CLIBS             := 
CXX               := $CXX
CXXFLAGS_COMPILE  := \$(FAST) -fmessage-length=0 -D_REENTRANT -D_THREAD_SAFE -Wall -Wno-char-subscripts
CXXLDFLAGS        := 
CXXLIBS           := 
CXXSHARED         := -Wl,-r -dynamic
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  FreeBSD-amd64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized
CC                := $CC
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -O3 -fPIC -pthread -D_REENTRANT -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -pthread -lthr
CXX               := $CXX
CXXFLAGS_COMPILE  := -O3 -fPIC -pthread -D_REENTRANT -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -pthread -lthr
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
CC                := $CC
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -g -pthread -D_REENTRANT -fPIC -Wall -Wno-char-subscripts -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs  
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -pthread -lthr
CXX               := $CXX
CXXFLAGS_COMPILE  := -g -pthread -D_REENTRANT -fPIC -Wall -Wno-char-subscripts -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -pthread -lthr
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  FreeBSD-amd64-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, debug, warnings
CC                := $CC
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -pg -O3 -pthread -D_REENTRANT -fPIC -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions
CLDFLAGS          := -pg -L/usr/local/lib
CLIBS             := -pthread -lthr
CXX               := $CXX
CXXFLAGS_COMPILE  := -pg -O3 -pthread -D_REENTRANT -fPIC -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions
CXXLDFLAGS        := -pg -L/usr/local/lib
CXXLIBS           := -pthread -lthr
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  Linux-i686)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux, optimized
CC                := $CC
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -O3 -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -pthread -ldl
CXX               := $CXX
CXXFLAGS_COMPILE  := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -O3 -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -pthread -ldl
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
CC                := $CC
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -m64 -fPIC -D_REENTRANT -O3 -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -pthread -ldl
CXX               := $CXX
CXXFLAGS_COMPILE  := -m64 -fPIC -D_REENTRANT -O3 -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -pthread -ldl
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
CC                := $CC
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -m64 -fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -pthread -ldl
CXX               := $CXX
CXXFLAGS_COMPILE  := -m64 -fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -pthread -ldl
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  Linux-amd64-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux64, optimized
CC                := $CC
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -pg -m64 -fPIC -D_REENTRANT -O3 -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -pthread -ldl
CXX               := $CXX
CXXFLAGS_COMPILE  := -pg -m64 -fPIC -D_REENTRANT -O3 -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -pthread -ldl
CXXSHARED         := -shared
ARFLAGS           := ruvs
INSTALL/          := $target/
EOF
    ;;
  Linux-ia64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux64, optimized
CC                := $CC
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -m64 -fPIC -D_REENTRANT -O3 -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -pthread -ldl
CXX               := $CXX
CXXFLAGS_COMPILE  := -m64 -fPIC -D_REENTRANT -O3 -D_THREAD_SAFE -pthread -fmessage-length=0 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -pthread -ldl
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
CC                := $CC
SHLIB_FLAGS       := -G #untested
CFLAGS_COMPILE    := -D_REENTRANT -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -Wall -Wno-char-subscripts -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := 
CLIBS             := -lpthread -lrt
CXX               := $CXX
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

echo "Configured."

#cat Make.compilers
