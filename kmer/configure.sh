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
CC              := cc
CFLAGS_ARCH     := -O3 -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS_ARCH   := -L/usr/local/lib
CLIBS_ARCH      := -llthread -llgcc_r
CXX             := g++
CXXFLAGS_ARCH   := -O3 -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS_ARCH := -L/usr/local/lib
CXXLIBS_ARCH    := -llthread -llgcc_r
ARFLAGS         := ruvs
EOF
    ;;
  freebsd-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized, profiled (same as FreeBSD, optimized, but includes -pg)
CC              := cc
CFLAGS_ARCH     := -pg -O3 -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS_ARCH   := -L/usr/local/lib
CLIBS_ARCH      := -llthread -llgcc_r
CXX             := g++
CXXFLAGS_ARCH   := -pg -O3 -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS_ARCH := -L/usr/local/lib
CXXLIBS_ARCH    := -llthread -llgcc_r
ARFLAGS         := ruvs
EOF
    ;;
  freebsd-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, debug, warnings
CC              := cc
CFLAGS_ARCH     := -g -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -Wshadow -Wtraditional -Wid-clash-16 -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs  
CLDFLAGS_ARCH   := -L/usr/local/lib
CLIBS_ARCH      := -llthread -llgcc_r
CXX             := g++
CXXFLAGS_ARCH   := -g -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -Wshadow -Wtraditional -Wid-clash-16 -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs  
CXXLDFLAGS_ARCH := -L/usr/local/lib
CXXLIBS_ARCH    := -llthread -llgcc_r
ARFLAGS         := ruvs
EOF
    ;;
  aix)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, optimized
#  Tested, verified on 30jan03
CC              := xlc_r
CFLAGS_ARCH     := -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict -qcpluscmt
CLDFLAGS_ARCH   := 
CLIBS_ARCH      := 
CXX             := xlC_r
CXXFLAGS_ARCH   := -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict
CXXLDFLAGS_ARCH := 
CXXLIBS_ARCH    := 
ARFLAGS         := -X 64 ruv
EOF
    ;;
  aix-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, optimized, profiled
#  Untested; might not work with optimization
CC              := xlc_r
CFLAGS_ARCH     := -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict -qcpluscmt
CLDFLAGS_ARCH   := 
CLIBS_ARCH      := 
CXX             := xlC_r
CXXFLAGS_ARCH   := -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict
CXXLDFLAGS_ARCH := 
CXXLIBS_ARCH    := 
ARFLAGS         := -X 64 ruv
EOF
    ;;
  aix-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, debug
ccDBUG          := -qheapdebug -qcheck=all -DMEMORY_DEBUG
ccDBUG          := 
CC              := xlc_r
CFLAGS_ARCH     := -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath -qcpluscmt
CLDFLAGS_ARCH   := 
CLIBS_ARCH      := 
CXX             := xlC_r
CXXFLAGS_ARCH   := -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath
CXXLDFLAGS_ARCH := 
CXXLIBS_ARCH    := 
ARFLAGS         := -X 64 ruv
EOF
    ;;
  compaq|tru64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Tru64, native compilers, optimized
CC              := cc
CFLAGS_ARCH     := -D_REENTRANT -pthread -w0 -fast
CLDFLAGS_ARCH   := 
CLIBS_ARCH      := -lpthread -lrt
CXX             := cxx
CXXFLAGS_ARCH   := -D_REENTRANT -pthread -w0 -fast
CXXLDFLAGS_ARCH := 
CXXLIBS_ARCH    := -lpthread -lrt
ARFLAGS         := ruv
CXX_TMP_ARCH    := cxx_repository
EOF
    ;;
  compaq-debug|tru64-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Tru64, native compilers, debug
CC              := cc
CFLAGS_ARCH     := -D_REENTRANT -pthread -g -w0 -trapuv
CLDFLAGS_ARCH   := 
CLIBS_ARCH      := -lpthread -lrt
CXX             := cxx
CXXFLAGS_ARCH   := -D_REENTRANT -pthread -gall -w0 -trapuv
CXXLDFLAGS_ARCH := 
CXXLIBS_ARCH    := -lpthread -lrt
ARFLAGS         := ruv
CXX_TMP_ARCH    := cxx_repository
EOF
    ;;
  linux)
    echo "Linux is UNTESTED!"
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux, optimized
CC              := cc
CFLAGS_ARCH     := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -pthread -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXX             := g++
CXXFLAGS_ARCH   := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -pthread -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
ARFLAGS         := ruv
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
CC              := gcc
CFLAGS_ARCH     := -D_REENTRANT -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS_ARCH   := 
CLIBS_ARCH      := -lpthread -lrt
CXX             := g++
CXXFLAGS_ARCH   := -D_REENTRANT -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer 
CXXLDFLAGS_ARCH := 
CXXLIBS_ARCH    := -lpthread -lrt
ARFLAGS         := ruv
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
    exit
    ;;
esac

cat <<EOF >> Make.compilers

CCDEP		:= gcc -MM -MG
CXXDEP		:= g++ -MM -MG

CFLAGS          :=
CXXFLAGS        :=
CLIBS		+= -lm
CXXLIBS         += -lm

EOF
