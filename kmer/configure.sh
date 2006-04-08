#!/bin/sh

#  Set up the build system -- need some symlinks to the build
#  directory.
#
if [ ! -e Makefile ] ; then
  if [ -e ../build/Makefile ] ; then
    ln -s ../build/Make.path .
    ln -s ../build/Make.rules .
    ln -s ../build/Make_utils .
    ln -s ../build/Makefile .
  elif [ -e build/Makefile ] ; then
    ln -s build/Make.path .
    ln -s build/Make.rules .
    ln -s build/Make_utils .
    ln -s build/Makefile .
  else
    echo "Hey, couldn't find the Makefile"
    exit
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
      target="osx$opts"
      ;;
    FreeBSD)
      target="freebsd$opts"
      ;;
    AIX)
      target="aix$opts"
      ;;
    OSF1)
      target="tru64$opts"
      ;;
    Linux)
      target="linux$opts"
      if [ `uname -m` = "x86_64" ] ; then
        target="linux64$opts"
      fi
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
#
#  Using this usually generates internal compiler errors: -mpowerpc64
#    FIXED with the "November 2004" update!
#  Using this breaks dynamic library building:  -mdynamic-no-pic
#    We could have instead included -fPIC on the compile line.
#
# -malign-natural changes the size of structures compared to -g.  -Wpadded supposedly warns when this happens
#
FAST              := -Wpadded -O3 -funroll-loops -fstrict-aliasing -fsched-interblock -falign-loops=16 -falign-jumps=16 -falign-functions=16 -falign-jumps-max-skip=15 -falign-loops-max-skip=15 -malign-natural -ffast-math -mpowerpc-gpopt -force_cpusubtype_ALL -fstrict-aliasing -mtune=G5 -mcpu=G5
FAST              := -fast -fPIC
CC                := gcc
SHLIB_FLAGS       := -dynamiclib
CFLAGS_COMPILE    := \$(FAST) -fmessage-length=0 -D_THREAD_SAFE -Wall
CLDFLAGS          := 
CLIBS             := 
CXX               := g++
CXXFLAGS_COMPILE  := \$(FAST) -fmessage-length=0 -D_THREAD_SAFE -Wall
CXXLDFLAGS        := 
CXXLIBS           := 
ARFLAGS           := ruvs
INSTALL/          := darwin/
EOF
    ;;
  osx-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  OS-X, debug
#
CC                := gcc
SHLIB_FLAGS       := -dynamiclib
CFLAGS_COMPILE    := -g3 -fmessage-length=0 -D_THREAD_SAFE -Wall
CLDFLAGS          := 
CLIBS             := 
CXX               := g++
CXXFLAGS_COMPILE  := -g3 -fmessage-length=0 -D_THREAD_SAFE -Wall
CXXLDFLAGS        := 
CXXLIBS           := 
ARFLAGS           := ruvs
INSTALL/          := darwin/
EOF
    ;;
  osx-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  OS-X, optimized, profiled (same as OS-X, optimized, but includes -pg)
#
FAST              := -pg -O3 -funroll-loops -fstrict-aliasing -fsched-interblock -falign-loops=16 -falign-jumps=16 -falign-functions=16 -falign-jumps-max-skip=15 -falign-loops-max-skip=15 -malign-natural -ffast-math -mpowerpc-gpopt -force_cpusubtype_ALL -fstrict-aliasing -mtune=G5 -mcpu=G5
CC                := gcc
SHLIB_FLAGS       := -dynamiclib
CFLAGS_COMPILE    := \$(FAST) -fmessage-length=0 -D_THREAD_SAFE -Wall
CLDFLAGS          := -pg
CLIBS             := 
CXX               := g++
CXXFLAGS_COMPILE  := \$(FAST) -fmessage-length=0 -D_THREAD_SAFE -Wall
CXXLDFLAGS        := -pg
CXXLIBS           := 
ARFLAGS           := ruvs
INSTALL/          := darwin/
EOF
    ;;
  freebsd)
    rm -f Make.compilers
    echo "Using linuxthreads by default!"
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized
THREADS           := -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads 
THREADL           := -llthread -llgcc_r
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -O3 \$(THREADS) -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL)
CXX               := g++
CXXFLAGS_COMPILE  := -O3 \$(THREADS) -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL)
ARFLAGS           := ruvs
INSTALL/          := freebsd/
EOF
    ;;
  freebsd-profile)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, optimized, profiled (same as FreeBSD, optimized, but includes -pg)
THREADS           := -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads 
THREADL           := -llthread -llgcc_r
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -pg -O3 \$(THREADS) -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL)
CXX               := g++
CXXFLAGS_COMPILE  := -pg -O3 \$(THREADS) -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL)
ARFLAGS           := ruvs
INSTALL/          := freebsd/
EOF
    ;;
  freebsd-debug)
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
CFLAGS_COMPILE    := -g \$(THREADS) -Wall -Wshadow -Wtraditional -Wpointer-arith -Wcast-qual \
  -Wcast-align -Wwrite-strings -Wconversion -Wstrict-prototypes -Wmissing-prototypes \
  -Wmissing-declarations -Wnested-externs  
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL)
CXX               := g++
CXXFLAGS_COMPILE  := -g \$(THREADS) -Wall -Wshadow -Wtraditional -Wpointer-arith -Wcast-qual \
  -Wcast-align -Wwrite-strings -Wconversion -Wstrict-prototypes -Wmissing-prototypes \
  -Wmissing-declarations  -Wnested-externs  
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL)
ARFLAGS           := ruvs
INSTALL/          := freebsd/
EOF
    ;;
  freebsd-debug-34)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  FreeBSD, debug, warnings
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -g \$(THREADS) -Wall -Wshadow -Wtraditional -Wpointer-arith -Wcast-qual -Wcast-align \
  -Wwrite-strings -Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL)
CXX               := g++
CXXFLAGS_COMPILE  := -g \$(THREADS) -Wall -Wshadow -Wtraditional -Wpointer-arith -Wcast-qual -Wcast-align \
  -Wwrite-strings -Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs  
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL)
ARFLAGS           := ruvs
INSTALL/          := freebsd/
EOF
    ;;
  aix)
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
ARFLAGS           := -X 64 ruv
INSTALL/          := aix/
EOF
    ;;
  aix-profile)
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
ARFLAGS           := -X 64 ruv
INSTALL/          := aix/
EOF
    ;;
  aix-debug)
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
ARFLAGS           := -X 64 ruv
INSTALL/          := aix/
EOF
    ;;
  aix-debug-full|aix-full-debug)
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
ARFLAGS           := -X 64 ruv
INSTALL/          := aix/
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
ARFLAGS           := ruv
CXX_TMP_ARCH      := cxx_repository
INSTALL/          := osf1/
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
ARFLAGS           := ruv
CXX_TMP_ARCH      := cxx_repository
INSTALL/          := osf1/
EOF
    ;;
  linux)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux, optimized
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -O3 \$(THREADS) -fmessage-length=0 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL) -ldl
CXX               := g++
CXXFLAGS_COMPILE  := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -O3 \$(THREADS) -fmessage-length=0 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL) -ldl
ARFLAGS           := ruvs
INSTALL/          := linux/
EOF
    ;;
  linux64)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux64, optimized
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -m64 -fPIC -D_REENTRANT -O3 \$(THREADS) -fmessage-length=0 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL) -ldl
CXX               := g++
CXXFLAGS_COMPILE  := -m64 -fPIC -D_REENTRANT -O3 \$(THREADS) -fmessage-length=0 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL) -ldl
ARFLAGS           := ruvs
INSTALL/          := linux64/
EOF
    ;;
  linux64-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux64, optimized
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -m64 -fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g \$(THREADS) -fmessage-length=0 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL) -ldl
CXX               := g++
CXXFLAGS_COMPILE  := -m64 -fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g \$(THREADS) -fmessage-length=0 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL) -ldl
ARFLAGS           := ruvs
INSTALL/          := linux64/
EOF
    ;;
  linux-debug)
    rm -f Make.compilers
    cat <<EOF > Make.compilers
# -*- makefile -*-
#  Linux, optimized
THREADS           := -D_THREAD_SAFE -pthread
THREADL           := -pthread
CC                := cc
SHLIB_FLAGS       := -shared
CFLAGS_COMPILE    := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g \$(THREADS) -fmessage-length=0 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := -L/usr/local/lib
CLIBS             := \$(THREADL) -ldl
CXX               := g++
CXXFLAGS_COMPILE  := -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_REENTRANT -g \$(THREADS) -fmessage-length=0 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CXXLDFLAGS        := -L/usr/local/lib
CXXLIBS           := \$(THREADL) -ldl
ARFLAGS           := ruvs
INSTALL/          := linux/
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
CFLAGS_COMPILE    := -D_REENTRANT -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
CLDFLAGS          := 
CLIBS             := -lpthread -lrt
CXX               := g++
CXXFLAGS_COMPILE  := -D_REENTRANT -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -Wall -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer 
CXXLDFLAGS        := 
CXXLIBS           := -lpthread -lrt
ARFLAGS           := ruv
INSTALL/          := solaris/
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
CLIBS             += -lm
CXXLIBS           += -lm
EOF

cat Make.compilers
