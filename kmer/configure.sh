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
    Darwin)
      target="osx"
      ;;
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
  osx)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  OS-X, optimized
CC                := gcc
CFLAGS_COMPILE    := -O3 -D_THREAD_SAFE -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := 
CLIBS             := 
CXX               := g++
CXXFLAGS_COMPILE  := -O3 -D_THREAD_SAFE -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := 
CXXLIBS           := 
ARFLAGS           := ruvs
EOF
    ;;
  freebsd)
    rm -f Make.compilers
    echo "Using linuxthreads by default!"
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized
CC                := cc
CFLAGS_COMPILE    := -O3 -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -llthread -llgcc_r
CXX               := g++
CXXFLAGS_COMPILE  := -O3 -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -llthread -llgcc_r
ARFLAGS           := ruvs
EOF
    ;;
  freebsd-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized, profiled (same as FreeBSD, optimized, but includes -pg)
CC                := cc
CFLAGS_COMPILE    := -pg -O3 -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -llthread -llgcc_r
CXX               := g++
CXXFLAGS_COMPILE  := -pg -O3 -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -llthread -llgcc_r
ARFLAGS           := ruvs
EOF
    ;;
  freebsd-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, debug, warnings
CC                := cc
CFLAGS_COMPILE    := -g -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -Wshadow -Wtraditional -Wid-clash-16 -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs  
CLDFLAGS          := -L/usr/local/lib
CLIBS             := -llthread -llgcc_r
CXX               := g++
CXXFLAGS_COMPILE  := -g -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads -Wall -Wshadow -Wtraditional -Wid-clash-16 -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs  
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := -llthread -llgcc_r
ARFLAGS           := ruvs
EOF
    ;;
  aix)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, optimized
#  Tested, verified on 30jan03
CC                := xlc_r
CFLAGS_COMPILE    := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict -qcpluscmt
CLDFLAGS          := 
CLIBS             := 
CXX               := xlC_r
CXXFLAGS_COMPILE  := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict
CXXLDFLAGS        := 
CXXLIBS           := 
ARFLAGS           := -X 64 ruv
EOF
    ;;
  aix-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, optimized, profiled
#  Untested; might not work with optimization
CC                := xlc_r
CFLAGS_COMPILE    := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict -qcpluscmt
CLDFLAGS          := 
CLIBS             := 
CXX               := xlC_r
CXXFLAGS_COMPILE  := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES -q64 -O3 -pg -qfullpath -qmaxmem=-1 -qarch=auto -qtune=auto -qcache=auto -qstrict
CXXLDFLAGS        := 
CXXLIBS           := 
ARFLAGS           := -X 64 ruv
EOF
    ;;
  aix-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, debug
ccDBUG            := 
CC                := xlc_r
CFLAGS_COMPILE    := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath -qcpluscmt
CLDFLAGS          := 
CLIBS             := 
CXX               := xlC_r
CXXFLAGS_COMPILE  := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath
CXXLDFLAGS        := 
CXXLIBS           := 
ARFLAGS           := -X 64 ruv
EOF
    ;;
  aix-debug-full|aix-full-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  AIX5, full debug, including memory checks
ccDBUG            := -qheapdebug -qinfo=uni:use -qcheck=all -qinitauto=FF -qflttrap=inv:ov:zero:en -qfloat=nans -DMEMORY_DEBUG
CC                := xlc_r
CFLAGS_COMPILE    := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath -qcpluscmt
CLDFLAGS          := \$(ccDBUG)
CLIBS             := 
CXX               := xlC_r
CXXFLAGS_COMPILE  := -qstaticinline -qthreaded -D_THREAD_SAFE -D_LARGE_FILES \$(ccDBUG) -q64 -g -qfullpath
CXXLDFLAGS        := \$(ccDBUG)
CXXLIBS           := 
ARFLAGS           := -X 64 ruv
EOF
    ;;
  compaq|tru64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Tru64, native compilers, optimized
CC                := cc
CFLAGS_COMPILE    := -D_REENTRANT -pthread -w0 -fast
CLDFLAGS          := 
CLIBS             := -lpthread -lrt
CXX               := cxx
CXXFLAGS_COMPILE  := -D_REENTRANT -pthread -w0 -fast
CXXLDFLAGS        := 
CXXLIBS           := -lpthread -lrt
ARFLAGS           := ruv
CXX_TMP_ARCH      := cxx_repository
EOF
    ;;
  compaq-debug|tru64-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Tru64, native compilers, debug
CC                := cc
CFLAGS_COMPILE    := -D_REENTRANT -pthread -g -w0 -trapuv
CLDFLAGS          := 
CLIBS             := -lpthread -lrt
CXX               := cxx
CXXFLAGS_COMPILE  := -D_REENTRANT -pthread -gall -w0 -trapuv
CXXLDFLAGS        := 
CXXLIBS           := -lpthread -lrt
ARFLAGS           := ruv
CXX_TMP_ARCH      := cxx_repository
EOF
    ;;
  linux)
    echo "Linux is UNTESTED!"
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux, optimized
CC                := cc
CFLAGS_COMPILE    := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -pthread -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXX               := g++
CXXFLAGS_COMPILE  := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -pthread -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
ARFLAGS           := ruv
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
CFLAGS_COMPILE    := -D_REENTRANT -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := 
CLIBS             := -lpthread -lrt
CXX               := g++
CXXFLAGS_COMPILE  := -D_REENTRANT -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer 
CXXLDFLAGS        := 
CXXLIBS           := -lpthread -lrt
ARFLAGS           := ruv
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
CCDEP		  := gcc -MM -MG
CXXDEP	          := g++ -MM -MG
CLIBS             += -lm
CXXLIBS           += -lm
EOF

cat Make.compilers
