#!/bin/sh

#  Script to generate Make.compilers, the configuration of the build.

#  XXX:  Add in efence support.  Need to add it into external

target=$1

#  If no target, try to figure out one based on uname.  This defaults to
#  the optimized target below.  If it works well, we can always use this
#  mechanism, and extend with "debug" or "profile" (e.g., "./configure.sh debug")
#
if [ "x$target" = "x" ] ; then
  case `uname` in
    FreeBSD)
      target="freebsd"
      ;;
    AIX)
      target="aix"
      ;;
    OSF1)
      target="tru64"
      ;;
    *)
      echo "ERROR: Unknown uname of `uname` -- try manual configuration."
      exit 1
      ;;
  esac
fi


case $target in
  freebsd)
    rm -f Make.compilers
    echo "Using linuxthreads by default!"
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized
THRINC  := -D_REENTRANT -pthread
THRLIB  := 
THRINC  := -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads
THRLIB  := -L/usr/local/lib -llthread -llgcc_r
CC      := cc
CFLAGS_ARCH  := -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CWARNS := 
CXX      := g++
CXXFLAGS_ARCH  := -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXWARNS := 
ARFLAGS  := ruvs
EOF
    ;;
  freebsd-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized, profiled
THRINC  := -D_REENTRANT -pthread
THRLIB  := 
CC      := cc
CFLAGS_ARCH  := -O3 -pg -Wall -funroll-loops -fexpensive-optimizations -finline-functions
CWARNS := 
CXX      := g++
CXXFLAGS_ARCH  := -O3 -pg -Wall -funroll-loops -fexpensive-optimizations -finline-functions
CXXWARNS := 
ARFLAGS  := ruvs
EOF
    ;;
  freebsd-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, debug, warnings
THRINC  := -D_REENTRANT -pthread
THRLIB  := 
CC      := cc
CFLAGS_ARCH  := -g -Wall
CWARNS := -Wshadow -Wtraditional -Wid-clash-16 -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs  
CXX      := g++
CXXFLAGS_ARCH  := -g -Wall
CXXWARNS := -Wshadow -Wtraditional -Wid-clash-16 -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs  
ARFLAGS  := ruvs
EOF
    ;;
  aix)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, optimized
#  Tested, verified on 30jan03
THRINC  := -qthreaded -D_THREAD_SAFE
THRLIB  := -lpthread
IBMLF	:= -D_LARGE_FILES
CC      := xlc_r
CFLAGS_ARCH  :=-q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict \$(IBMLF) -qcpluscmt
CWARNS := 
CXX      := xlC_r
CXXFLAGS_ARCH  := -q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict \$(IBMLF)
CXXWARNS :=
ARFLAGS  := -X 64 ruv
EOF
    ;;
  aix-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, optimized, profiled
#  Untested; might not work with optimization
THRINC  := -qthreaded -D_THREAD_SAFE
THRLIB  := -lpthread
IBMLF	:= -D_LARGE_FILES
CC      := xlc_r
CFLAGS_ARCH  :=  -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict \$(IBMLF) -qcpluscmt
CWARNS := 
CXX      := xlC_r
CXXFLAGS_ARCH  := -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict \$(IBMLF)
CXXWARNS :=
ARFLAGS  := -X 64 ruv
EOF
    ;;
  aix-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, debug
THRINC  := -qthreaded -D_THREAD_SAFE
THRLIB  := -lpthread
IBMLF	:= -D_LARGE_FILES
CC      := xlc_r
CFLAGS_ARCH  := -q64 -g -qfullpath \$(IBMLF) -qcpluscmt
CWARNS := 
CXX      := xlC_r
CXXFLAGS_ARCH  := -q64 -g -qfullpath \$(IBMLF)
CXXWARNS :=
ARFLAGS  := -X 64 ruv
EOF
    ;;
  compaq|tru64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Tru64, native compilers, optimized
THRINC  := 
THRLIB  := -lpthread -lrt
CC      := cc
CFLAGS_ARCH  := -D_REENTRANT -pthread -w0 -fast
CWARNS := 
CXX      := cxx
CXXFLAGS_ARCH  := -D_REENTRANT -pthread -w0 -fast
CXXWARNS :=
ARFLAGS  := ruv

CXX_TMP_ARCH := cxx_repository

EOF
    ;;
  compaq-debug|tru64-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Tru64, native compilers, debug
THRINC  := -D_REENTRANT
THRLIB  := -pthread -lpthread -lrt
CC      := cc
CFLAGS_ARCH  := -D_REENTRANT -pthread -g -w0 -trapuv
CWARNS := 
CXX      := cxx
CXXFLAGS_ARCH  := -D_REENTRANT -pthread -gall -w0 -trapuv
CXXWARNS :=
ARFLAGS  := ruv

CXX_TMP_ARCH := cxx_repository

EOF
    ;;
  linux)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux, optimized
THRINC  := -D_REENTRANT -pthread
THRLIB  := 
CC      := cc
CFLAGS_ARCH  := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CWARNS := 
CXX      := g++
CXXFLAGS_ARCH  := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXWARNS := 
ARFLAGS  := ruv
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
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Solaris, gcc optimized
#
THRINC  := -D_REENTRANT
THRLIB  := -lpthread -lrt
SUNLF	:= -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CC      := gcc
CFLAGS_ARCH  := -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer \$(SUNLF)
CWARNS := 
CXX      := g++
CXXFLAGS_ARCH  := -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer \$(SUNLF)
CXXWARNS := 
ARFLAGS  := ruv
EOF
    ;;
  solaris-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Solaris, gcc debug
THRINC  := -D_REENTRANT
THRLIB  := -lpthread -lrt
SUNLF	:= -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CC      := gcc
CFLAGS_ARCH  := -g -Wall \$(SUNLF)
CWARNS := 
CXX      := g++
CXXFLAGS_ARCH  := -g -Wall \$(SUNLF)
CXXWARNS := 
ARFLAGS  := ruv
EOF
    ;;
  *)
    echo "usage: $0 <configuration>"
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
    echo "          linux                     Linux, optimized (STALE)"
    echo ""
    echo "          solaris                   Solaris, gcc, optimized (STALE)"
    echo "          solaris-debug             Solaris, gcc, debug and warnings (STALE)"
    exit
    ;;
esac

cat <<EOF >> Make.compilers

CFLAGS_ARCH   += ${CWARNS}
CXXFLAGS_ARCH += ${CXXWARNS}

CCDEP		 :=gcc -MM -MG
CXXDEP		 :=g++ -MM -MG

CFLAGS          :=
CXXFLAGS        :=
CLIBS		+=-lm
CXXLIBS         +=-lm

EOF