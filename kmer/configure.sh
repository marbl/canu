#!/bin/sh

#  Script to generate Makefile-configuration, the configuration of the build.

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
    rm -f Makefile-configuration
    echo "Using linuxthreads by default!"
    cat <<EOF > Makefile-configuration
#  FreeBSD, optimized
THRINC  = -D_REENTRANT -pthread
THRLIB  = 
THRINC  = -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads
THRLIB  = -L/usr/local/lib -llthread -llgcc_r
cc      = cc
ccOPTS  = -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
ccWARNS = 
CC      = g++
CCOPTS  = -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CCWARNS = 
AROPTS  = ruvs
EOF
    ;;
  freebsd-profile)
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  FreeBSD, optimized, profiled
THRINC  = -D_REENTRANT -pthread
THRLIB  = 
cc      = cc
ccOPTS  = -O3 -pg -Wall -funroll-loops -fexpensive-optimizations -finline-functions
ccWARNS = 
CC      = g++
CCOPTS  = -O3 -pg -Wall -funroll-loops -fexpensive-optimizations -finline-functions
CCWARNS = 
AROPTS  = ruvs
EOF
    ;;
  freebsd-debug)
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  FreeBSD, debug, warnings
THRINC  = -D_REENTRANT -pthread
THRLIB  = 
cc      = cc
ccOPTS  = -g -Wall
ccWARNS = -Wshadow -Wtraditional -Wid-clash-16 -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs  
CC      = g++
CCOPTS  = -g -Wall
CCWARNS = -Wshadow -Wtraditional -Wid-clash-16 -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs  
AROPTS  = ruvs
EOF
    ;;
  aix)
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  AIX5, optimized
#  Tested, verified on 30jan03
THRINC  = -qthreaded -D_THREAD_SAFE
THRLIB  = -lpthread
IBMLF	= -D_LARGE_FILES
cc      = xlc
ccOPTS  = -q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict \$(IBMLF) -qcpluscmt
ccWARNS = 
CC      = xlC
CCOPTS  = -q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict \$(IBMLF)
CCWARNS =
AROPTS  = -X 64 ruv
EOF
    ;;
  aix-profile)
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  AIX5, optimized, profiled
#  Untested; might not work with optimization
THRINC  = -qthreaded -D_THREAD_SAFE
THRLIB  = -lpthread
IBMLF	= -D_LARGE_FILES
cc      = xlc
ccOPTS  = -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict \$(IBMLF) -qcpluscmt
ccWARNS = 
CC      = xlC
CCOPTS  = -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict \$(IBMLF)
CCWARNS =
AROPTS  = -X 64 ruv
EOF
    ;;
  aix-debug)
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  AIX5, debug
THRINC  = -qthreaded -D_THREAD_SAFE
THRLIB  = -lpthread
IBMLF	= -D_LARGE_FILES
cc      = xlc
ccOPTS  = -q64 -g -qfullpath \$(IBMLF) -qcpluscmt
ccWARNS = 
CC      = xlC
CCOPTS  = -q64 -g -qfullpath \$(IBMLF)
CCWARNS =
AROPTS  = -X 64 ruv
EOF
    ;;
  compaq|tru64)
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  Tru64, native compilers, optimized
THRINC  = -D_REENTRANT
THRLIB  = -pthread -lpthread -lrt
cc      = cc
ccOPTS  = -w0 -fast
ccWARNS = 
CC      = cxx
CCOPTS  = -w0 -fast
CCWARNS =
AROPTS  = ruv
EOF
    ;;
  compaq-debug|tru64-debug)
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  Tru64, native compilers, debug
THRINC  = -D_REENTRANT
THRLIB  = -pthread -lpthread -lrt
cc      = cc
ccOPTS  = -g -w0 -trapuv
ccWARNS = 
CC      = cxx
CCOPTS  = -gall -w0 -trapuv
CCWARNS =
AROPTS  = ruv
EOF
    ;;
  linux)
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  Linux, optimized
THRINC  = -D_REENTRANT -pthread
THRLIB  = 
cc      = cc
ccOPTS  = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
ccWARNS = 
CC      = g++
CCOPTS  = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CCWARNS = 
AROPTS  = ruv
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
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  Solaris, gcc optimized
#
THRINC  = -D_REENTRANT
THRLIB  = -lpthread -lrt
SUNLF	= -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
cc      = gcc
ccOPTS  = -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer \$(SUNLF)
ccWARNS = 
CC      = g++
CCOPTS  = -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer \$(SUNLF)
CCWARNS = 
AROPTS  = ruv
EOF
    ;;
  solaris-debug)
    rm -f Makefile-configuration
    cat <<EOF > Makefile-configuration
#  Solaris, gcc debug
THRINC  = -D_REENTRANT
THRLIB  = -lpthread -lrt
SUNLF	= -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
cc      = gcc
ccOPTS  = -g -Wall \$(SUNLF)
ccWARNS = 
CC      = g++
CCOPTS  = -g -Wall \$(SUNLF)
CCWARNS = 
AROPTS  = ruv
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

cat <<EOF >> Makefile-configuration

.c.o: \$*.c
	\$(cc) \$(ccOPTS) \$(ccWARNS) \$(INCLUDE) -c -o \$*.o \$*.c
.C.o: \$*.C
	\$(CC) \$(CCOPTS) \$(CCWARNS) \$(INCLUDE) -c -o \$*.o \$*.C
EOF
