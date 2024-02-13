# -*- Makefile -*-
#
# boilermake: A reusable, but flexible, boilerplate Makefile.
#
# Copyright 2008, 2009, 2010 Dan Moulding, Alan T. DeKok
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#
# Significantly modified for Marbl projects.  Get the original from
#   https://github.com/dmoulding/boilermake
#

include Makefile.boilermake

#
#  Check that Make and git are compatible, fail if not.
#

define TEST_VERSION
  path   := $$(shell which ${1})
  vers   := $$(shell ${1}     --version | head -n 1 | cut -d\  -f ${2})
  versn  := $$(shell echo $${vers} | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$$$/&00/')
  min    := $$(shell echo ${3}     | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$$$/&00/')

  ifeq ($$(shell expr $${versn} \>= $${min}), 0)
    $$(error '$${path}' version '$${vers}' too old; at least version $3 is required)
  endif
endef

$(eval $(call TEST_VERSION, git,     3, 2.12))
$(eval $(call TEST_VERSION, ${MAKE}, 3, 3.81))

#
#  Initialize submodules if they aren't here yet.
#   - if the submodule directory exists and it is empty, then have
#     git fetch the submodule (with some fancy logging).
#

define INIT_SUBMODULE
  ifeq ($$(wildcard ${1}),${1})
  ifeq ($$(wildcard ${1}/*),)
    $$(info Fetching submodule '${1}')
    $$(shell git submodule update --init ${2} ${1} 2>&1 | awk '{ print " - "$$$$0 }' 1>&2)
    $$(info )
  endif
  endif
endef

$(eval $(call INIT_SUBMODULE,utility,))
$(eval $(call INIT_SUBMODULE,meryl,))
$(eval $(call INIT_SUBMODULE,seqrequester,))

#
#  Set compiler and flags based on operating system and compiler toolchain.
#  This must be called AFTER the compiler is chosen, and after the decision
#  to include stack trace support is make.
#
#  User-supplied flags are handled:
#    If a parameter to make ("make CXXFLAGS=whatever"):
#      flags are used exclusively; no default flags are supplied
#    If an environment variable ("CXXFLAGS=whatever make"):
#      flags are appended to the end of the default flags.
#
#  Parameters to the function:
#    $1 - extra flags for CFLAGS/CXXFLAGS (e.g., -I/path)
#    $2 - extra flags for LDFLAGS (e.g., -L/path)
#
#  Options:
#    BUILDOPTIMIZED  will disable debug symbols (leaving it just optimized)
#    BUILDDEBUG      will disable optimization  (leaving it just with debug symbols)
#    BUILDSTACKTRACE will enable stack trace on crashes, only works for Linux
#                    set to 0 on command line to disable (it's enabled by default for Linux)
#
#  BUILDPROFILE used to add -pg to LDFLAGS, and remove -D_GLIBCXX_PARALLE from CXXFLAGS and LDFLAGS,
#  and remove -fomit-frame-pointer from CXXFLAGS.  It added a bunch of complication and wasn't
#  really used.
#
#  BUILDJEMALLOC will enable jemalloc library support.
#
define SET_DEFAULT_FLAGS
  CFLAGSUSER   := $${CFLAGS}
  CXXFLAGSUSER := $${CXXFLAGS}
  LDFLAGSUSER  := $${LDFLAGS}

  #  Figure out which tool chain we're using and set up initial flags.
  #   - GCC seems to want '-dumpfullversion' instead of '-dumpversion'
  #     to return the full X.Y.Z version number.
  #   - GCC has deprecated c++2a but gcc-9 doesn't know the correct
  #     option c++20.

  ifeq ($$(shell echo `$$(CXX) --version 2>&1 | grep -c clang`), 0)
    TOOLCHAIN  := GNU
    TOOLCHAINV := $$(shell ${CXX} -dumpfullversion)
    CFLAGS     := -pthread -fopenmp -fPIC -std=c17   ${1}
    CXXFLAGS   := -pthread -fopenmp -fPIC -std=c++2a ${1}
    LDFLAGS    := -pthread -fopenmp                  ${2} -lm
  else
    TOOLCHAIN  := LLVM
    TOOLCHAINV := $$(shell ${CXX} -dumpversion)
    CFLAGS     := -pthread -fopenmp -fPIC -std=c17   ${1}
    CXXFLAGS   := -pthread -fopenmp -fPIC -std=c++20 ${1}
    LDFLAGS    := -pthread -fopenmp                  ${2} -lm
  endif

  #  Build up our preferred set of flags.

  ifeq ($$(BUILDOPTIMIZED), 1)   #  Enable debugging and add debug-aware
  else                           #  optimization unless BUILDOPTIMIZED is set.
    CXXFLAGS += -Og -g3
    CFLAGS   += -Og -g3
  endif

  ifeq ($$(BUILDDEBUG), 1)       #  Add full optimization unless BUILDDEBUG is set.
  else ifeq ($$(TOOLCHAIN), GNU)
    CXXFLAGS += -O4 -funroll-loops -finline-functions -fomit-frame-pointer -fexpensive-optimizations
    CFLAGS   += -O4 -funroll-loops -finline-functions -fomit-frame-pointer -fexpensive-optimizations
  else ifeq ($$(TOOLCHAIN), LLVM)
    CXXFLAGS += -O3 -funroll-loops -finline-functions -fomit-frame-pointer
    CFLAGS   += -O3 -funroll-loops -finline-functions -fomit-frame-pointer
  else
    CXXFLAGS += -O3
    CFLAGS   += -O3
  endif

  ifeq ($${BUILDSTACKTRACE}, 1)  #  Add (or disable) stacktrace support.
    CXXFLAGS += -DLIBBACKTRACE
    CFLAGS   += -DLIBBACKTRACE
  else
    CXXFLAGS += -DNOBACKTRACE
    CFLAGS   += -DNOBACKTRACE
  endif

  ifeq ($$(BUILDJEMALLOC), 1)    #  Add jemalloc support.
    CXXFLAGS += -DJEMALLOC -I`jemalloc-config --includedir`
    LDFLAGS  += -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -ljemalloc `jemalloc-config --libs`
  endif

  ifeq ($$(MACHINETYPE), amd64)  #  Add SIMD support (htslib and meryl-utility need this)
    CFLAGS   += -mxsave
    CXXFLAGS += -mxsave
  endif

  #  Enable gobs of warnings, but then disable some of the annoying ones, a
  #  few of which are compiler dependent.

  DEFAULTWARNINGS  =
  #EFAULTWARNINGS += -DMACHINETYPE_${MACHINETYPE}
  #EFAULTWARNINGS += -pedantic
  DEFAULTWARNINGS += -Wall -Wextra -Wformat
  DEFAULTWARNINGS += -Wno-char-subscripts
  DEFAULTWARNINGS += -Wno-sign-compare
  DEFAULTWARNINGS += -Wno-unused-function
  DEFAULTWARNINGS += -Wno-unused-parameter
  DEFAULTWARNINGS += -Wno-unused-variable
  DEFAULTWARNINGS += -Wno-deprecated-declarations
  #EFAULTWARNINGS += -Wno-format -Wno-format-truncation  #  Disable 'sprintf() into possibly smaller buffer'

  ifeq ($$(TOOLCHAIN), LLVM)
    CFLAGS   += $${DEFAULTWARNINGS} -Wno-format
    CXXFLAGS += $${DEFAULTWARNINGS} -Wno-format
  endif

  ifeq ($$(TOOLCHAIN), GNU)
    CFLAGS   += $${DEFAULTWARNINGS} -Wno-format-truncation
    CXXFLAGS += $${DEFAULTWARNINGS} -Wno-format-truncation
  endif

  DEFAULTWARNINGS =

  #  Append user-supplied flags, if any exist.

  ifneq ($${CFLAGSUSER}, )
    CFLAGS   += $${CFLAGSUSER}
  endif
  ifneq ($${CXXFLAGSUSER}, )
    CXXFLAGS += $${CXXFLAGSUSER}
  endif
  ifneq ($${LDFLAGSUSER}, )
    LDFLAGS += $${LDFLAGSUSER}
  endif
endef


#
#  Linux config.
#
ifeq (${OSTYPE}, Linux)
  CC        ?= gcc
  CXX       ?= g++

  BUILDSTACKTRACE ?= 1

  $(eval $(call SET_DEFAULT_FLAGS,,))
endif


#
#  FreeBSD config.
#
#  If building in the FreeBSD ports system, use the architecture as defined
#  there (technically, -p, not -m) and assume compiler and most options are
#  already defined correctly.
#
#  We used to (for non-ports builds) default to gcc, mostly because the rpath
#  was necessary.  But clang is working well, and rpath doesn't seem to be
#  necessary on FreeBSD13.
#    ifeq ($(origin CXX), default)
#      CC    = gcc9
#      CXX   = g++9
#      CCLIB = -rpath /usr/local/lib/gcc9
#    endif
#
ifeq (${OSTYPE}, FreeBSD)
ifeq (${CANU_BUILD_ENV}, ports)
  MACHINETYPE=${ARCH}

  CXXFLAGS  += -pthread -fopenmp -fPIC
  CFLAGS    += -pthread -fopenmp -fPIC
  LDFLAGS   += -pthread -fopenmp -lm
else
  CC        ?= cc
  CXX       ?= c++

  BUILDSTACKTRACE ?= 1

  $(eval $(call SET_DEFAULT_FLAGS,-I/usr/local/include,-L/usr/local/lib -lexecinfo))
endif
endif


#
#  MacOS config.
#
#  The default compiler in MacOS _still_ doesn't support OpenMP, so
#  we try a bunch of common alternate compiler names and use the first
#  one that exists.
#
#  If from MacPorts:             If from homebrew:
#    port install gcc9             brew install gcc
#    port select gcc mp-gcc9       brew install llvm
#
#  Homebrew calls its binaries clang and clang++ and uses directories to
#  differentiate.  (While there is clang-16, there is no clang++-16, and we
#  use clang/clang++ for consistency.)
#    /opt/homebrew/opt/gcc@11/bin/gcc-12   /opt/homebrew/opt/llvm@16/bin/clang
#    /opt/homebrew/opt/gcc@11/bin/g++-12   /opt/homebrew/opt/llvm@16/bin/clang++
#
#  MacPorts puts its binaries in the global bin directory with a suffix to
#  differentiate.
#    /opt/local/bin/gcc-mp-16           /opt/local/bin/clang-mp-16
#    /opt/local/bin/g++-mp-16           /opt/local/bin/clang++-mp-16
#
define FIND_MACOS_COMPILER
  ifeq ($$(CC), cc)
  ifneq ($$(wildcard $${BREW}/opt/gcc@${1}/bin/gcc-${1}), )
  ifneq ($$(wildcard $${BREW}/opt/gcc@${1}/bin/g++-${1}), )
    #$$(info Detected gcc-${1} installed via Homebrew.)
    CC=$$(abspath $${BREW}/opt/gcc@${1}/bin/gcc-${1})
    CXX=$$(abspath $${BREW}/opt/gcc@${1}/bin/g++-${1})
  endif
  endif
  endif

  ifeq ($$(CC), cc)
  ifneq ($$(wildcard $${BREW}/opt/llvm@${1}/bin/clang), )
  ifneq ($$(wildcard $${BREW}/opt/llvm@${1}/bin/clang++), )
    #$$(info Detected llvm-${1} installed via Homebrew.)
    CC=$$(abspath $${BREW}/opt/llvm@${1}/bin/clang)
    CXX=$$(abspath $${BREW}/opt/llvm@${1}/bin/clang++)
  endif
  endif
  endif

  ifeq ($$(CC), cc)
  ifneq ($$(wildcard $${PORT}/bin/gcc-mp-${1}), )
  ifneq ($$(wildcard $${PORT}/bin/g++-mp-${1}), )
    #$$(info Detected gcc-${1} installed via MacPorts.)
    CC=$$(abspath $${PORT}/bin/gcc-mp-${1})
    CXX=$$(abspath $${PORT}/bin/g++-mp-${1})
  endif
  endif
  endif

  ifeq ($$(CC), cc)
  ifneq ($$(wildcard $${PORT}/bin/clang-mp-${1}), )
  ifneq ($$(wildcard $${PORT}/bin/clang++-mp-${1}), )
    #$$(info Detected llvm-${1} installed via MacPorts.)
    CC=$$(abspath $${PORT}/bin/clang-mp-${1})
    CXX=$$(abspath $${PORT}/bin/clang++-mp-${1})
  endif
  endif
  endif
endef

ifeq (${OSTYPE}, Darwin)
  BREW := $(abspath $(lastword $(shell brew 2>/dev/null config | grep HOMEBREW_PREFIX)))
  PORT := $(abspath $(dir $(abspath $(dir $(shell which port)))))

  #$(info Detected Homebrew in '${BREW}')
  #$(info Detected MacPorts in '${PORT}')

  $(eval $(call FIND_MACOS_COMPILER,20))
  $(eval $(call FIND_MACOS_COMPILER,19))
  $(eval $(call FIND_MACOS_COMPILER,18))
  $(eval $(call FIND_MACOS_COMPILER,17))
  $(eval $(call FIND_MACOS_COMPILER,16))
  $(eval $(call FIND_MACOS_COMPILER,15))
  $(eval $(call FIND_MACOS_COMPILER,14))
  $(eval $(call FIND_MACOS_COMPILER,13))
  $(eval $(call FIND_MACOS_COMPILER,12))
  $(eval $(call FIND_MACOS_COMPILER,11))
  $(eval $(call FIND_MACOS_COMPILER,10))
  $(eval $(call FIND_MACOS_COMPILER,9))
  $(eval $(call FIND_MACOS_COMPILER,8))

  $(eval $(call SET_DEFAULT_FLAGS,,))
endif


#
#  Cygwin config.
#
ifneq (,$(findstring CYGWIN, ${OSTYPE}))
  CC        ?= gcc
  CXX       ?= g++

  $(eval $(call SET_DEFAULT_FLAGS,,))
endif


#
#  Test a small compile.
#
#  This has only failed on MacOS with default compilers that do not support
#  -fopenmp.
#
#  Note that .SHELLSTATUS was introduced in Make 4.2 but stupid MacOS ships
#  with make 3.81, so we need to do the success/fail test the hard way.  On
#  the otherhand, this does mean we can clean up temporary files right after.
#

COMPILETEST := $(shell echo "int main(void) { return 0; }" | ${CXX} ${CXXFLAGS} ${LDFLAGS} -x c++ -o /tmp/test-fopenmp - && echo pass ; rm -f /tmp/test-fopenmp /tmp/test-fopenmp.C)

ifneq ($(COMPILETEST), pass)
  $(warning )
  $(warning Unable to compile OpenMP programs with with:)
  $(warning $_    ${TOOLCHAIN} ${CXX} ${TOOLCHAINV}.)
  $(warning $_  CXXFLAGS flags:)
  $(foreach FLAG,${CXXFLAGS},$(warning $_    ${FLAG}))
  $(warning $_  LDFLAGS flags:)
  $(foreach FLAG,${LDFLAGS},$(warning $_    ${FLAG}))
  $(warning Please install GCC or LLVM with OpenMP support and)
  $(warning specify this compiler on the command line, e.g.,)
  $(warning $_    make CC=/path/to/gcc CXX=/path/to/g++)
  $(warning )
  $(error Unsupported compiler)
endif

#
# Define the two exported targets, 'all' and 'clean' as phony rules that
# point to the real rules.  This is done before including main.mk so that it
# can define any weird one-off special-case targets without them becoming the
# first (and thus default) rule.
#
# Then recursively include user-supplied submakefiles.
#
# Then define the real 'all' and 'clean' targets.
#
# Targets 'doall' and 'doclean' are used internally to add more targets to
# the global build.  'doclean-final' is a trouble maker because it needs to
# run after all the other 'doclean' steps, and the only way to ensure that is
# to give it a new target.

.PHONY: all
all: doall

.PHONY: clean
clean: doclean-final

$(eval $(call INCLUDE_SUBMAKEFILE,main.mk))

.PHONY: doall
doall: $(addprefix ${TARGET_DIR}/,${ALL_TGTS})
	@echo ""
	@echo "Success!"
	@echo "${MODULE} installed in ${TARGET_DIR}/bin/${MODULE}"
	@echo ""

.PHONY: doclean-final
doclean-final: doclean
	if [ -d ${TARGET_DIR} ] ; then find ${TARGET_DIR} -type d -print | sort -r | xargs -n 100 rmdir ; fi
	@echo ""
	@echo "Cleaned."
	@echo ""

.PHONY: doclean
doclean:

#
#  Let boilermake do its thing.
#

$(eval $(call BOILERMAKE))

#
#  Generate a version number.  This needs to come AFTER submakefiles are
#  included (so we know what ${MODULE} is) but before we build anything
#  (because this generates 'version.H' and also prints information about
#  submodule versions).
#
#  But, because make helpfully squashes all the lines to a single line we
#  need to do something funky and write the version information for C++ and
#  for make to one file - which then is included here and in source files.
#  Makefile.boilermake adds this file as a dependency to ALL source files.
#
#  For projects that include meryl-utility as a submodule (the usual case)
#  write version.H in the submodule src/ directory; for meryl-utility itself,
#  write in the current directory.
#   - This _should_ be indentifying if the submodule exists in
#     utility/src, but on the very first build (when the submodules
#     are initialized above) a file-based detection was reliably failing,
#     as if `ifeq ($(wildcard utility/src), utility/src)` was tested
#     before the submodule was init'd.  Rather than deal with all
#     that baloney, we just assume that any build using this build
#     system is either meryl-utility itself or something that uses it.
#

#$(info Calling  ./${VERSION_H:.H=.pl} ${MODULE} ${VERSION} ${VERSION_H})
VERSION_H := $(shell ./${VERSION_H:.H=.pl} ${MODULE} ${VERSION} ${VERSION_H})

ifeq (${VERSION_H}, )
$(error Failed to create version.H, possibly missing the 'meryl-utility' submodule)
endif

include ${VERSION_H}

#
#  Log what compiler we're using and start the build.
#

$(info ${BUILDING})
$(info For '${OSTYPE}' '${OSVERSION}' as '${MACHINETYPE}' into '${TARGET_DIR}/{bin,obj}'.)
$(info Using ${TOOLCHAIN} '$(shell which ${CXX})' version '${TOOLCHAINV}'.)
ifneq (${CFLAGSUSER},)
$(info Using user-supplied CFLAGS   '${CFLAGSUSER}'.)
endif
ifneq (${CXXFLAGSUSER},)
$(info Using user-supplied CXXFLAGS '${CXXFLAGSUSER}'.)
endif
ifneq (${LDFLAGSUSER},)
$(info Using user-supplied LDFLAGS  '${LDFLAGSUSER}'.)
endif
$(info )
