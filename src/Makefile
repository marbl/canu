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

# Caution: Don't edit this Makefile! Create your own main.mk and other
#          submakefiles, which will be included by this Makefile.
#          Only edit this if you need to modify boilermake's behavior (fix
#          bugs, add features, etc).

# Note: Parameterized "functions" in this makefile that are marked with
#       "USE WITH EVAL" are only useful in conjuction with eval. This is
#       because those functions result in a block of Makefile syntax that must
#       be evaluated after expansion. Since they must be used with eval, most
#       instances of "$" within them need to be escaped with a second "$" to
#       accomodate the double expansion that occurs when eval is invoked.

#  Before doing ANYTHING, initialize submodules...if the version of
#  git is compatible.
gitv := $(shell git --version | cut -d\  -f 3 | cut -c 1)
ifeq (1, $(gitv))
  gitv := $(shell git --version | cut -d\  -f 3)
  $(error git '$(shell which git)' version '$(gitv)' too old; at least version 2.12 is required)
endif

ifeq ($(wildcard utility/src/Makefile), )
  $(info $(shell git submodule update --init utility))
  $(info $(space))
endif
ifeq ($(wildcard meryl/src/Makefile), )
  $(info $(shell git submodule update --init meryl))
  $(info $(space))
endif
ifeq ($(wildcard seqrequester/src/Makefile), )
  $(info $(shell git submodule update --init seqrequester))
  $(info $(space))
endif

# ADD_CLEAN_RULE - Parameterized "function" that adds a new rule and phony
#   target for cleaning the specified target (removing its build-generated
#   files).
#
#   USE WITH EVAL
#
define ADD_CLEAN_RULE
    clean: clean_${1}
    .PHONY: clean_${1}
    clean_${1}:
	$$(strip rm -f ${TARGET_DIR}/${1} $${${1}_OBJS:%.o=%.[doP]})
	$${${1}_POSTCLEAN}
endef

# ADD_OBJECT_RULE - Parameterized "function" that adds a pattern rule for
#   building object files from source files with the filename extension
#   specified in the second argument. The first argument must be the name of the
#   base directory where the object files should reside (such that the portion
#   of the path after the base directory will match the path to corresponding
#   source files). The third argument must contain the rules used to compile the
#   source files into object code form.
#
#   USE WITH EVAL
#
define ADD_OBJECT_RULE
${1}/%.o: ${2} utility/src/utility/version.H
	${3}
endef

# ADD_TARGET_RULE - Parameterized "function" that adds a new target to the
#   Makefile. The target may be an executable or a library. The two allowable
#   types of targets are distinguished based on the name: library targets must
#   end with the traditional ".a" extension.
#
#   USE WITH EVAL
#
define ADD_TARGET_RULE
    ifeq "$$(suffix ${1})" ".a"
        # Add a target for creating a static library.
        $${TARGET_DIR}/${1}: $${${1}_OBJS}
	        @mkdir -p $$(dir $$@)
	        $$(strip $${AR} $${ARFLAGS} $$@ $${${1}_OBJS})
	        $${${1}_POSTMAKE}
    else
      # Add a target for linking an executable. First, attempt to select the
      # appropriate front-end to use for linking. This might not choose the
      # right one (e.g. if linking with a C++ static library, but all other
      # sources are C sources), so the user makefile is allowed to specify a
      # linker to be used for each target.
      ifeq "$$(strip $${${1}_LINKER})" ""
          # No linker was explicitly specified to be used for this target. If
          # there are any C++ sources for this target, use the C++ compiler.
          # For all other targets, default to using the C compiler.
          ifneq "$$(strip $$(filter $${CXX_SRC_EXTS},$${${1}_SOURCES}))" ""
              ${1}_LINKER = $${CXX}
          else
              ${1}_LINKER = $${CC}
          endif
      endif

      $${TARGET_DIR}/${1}: $${${1}_OBJS} $${${1}_PREREQS}
	      @mkdir -p $$(dir $$@)
	      $$(strip $${${1}_LINKER} -o $$@ $${LDFLAGS} $${${1}_LDFLAGS} $${${1}_OBJS} $${${1}_LDLIBS} $${LDLIBS})
	      $${${1}_POSTMAKE}
  endif
endef

# CANONICAL_PATH - Given one or more paths, converts the paths to the canonical
#   form. The canonical form is the path, relative to the project's top-level
#   directory (the directory from which "make" is run), and without
#   any "./" or "../" sequences. For paths that are not  located below the
#   top-level directory, the canonical form is the absolute path (i.e. from
#   the root of the filesystem) also without "./" or "../" sequences.
define CANONICAL_PATH
$(patsubst ${CURDIR}/%,%,$(abspath ${1}))
endef

# COMPILE_C_CMDS - Commands for compiling C source code.
define COMPILE_C_CMDS
	@mkdir -p $(dir $@)
	$(strip ${CC} -o $@ -c -MD ${CFLAGS} ${SRC_CFLAGS} ${INCDIRS} \
	    ${SRC_INCDIRS} ${SRC_DEFS} ${DEFS} $<)
	@cp ${@:%$(suffix $@)=%.d} ${@:%$(suffix $@)=%.P}; \
	 sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	     -e '/^$$/ d' -e 's/$$/ :/' < ${@:%$(suffix $@)=%.d} \
	     >> ${@:%$(suffix $@)=%.P}; \
	 rm -f ${@:%$(suffix $@)=%.d}
endef

# COMPILE_CXX_CMDS - Commands for compiling C++ source code.
define COMPILE_CXX_CMDS
	@mkdir -p $(dir $@)
	$(strip ${CXX} -o $@ -c -MD ${CXXFLAGS} ${SRC_CXXFLAGS} ${INCDIRS} \
	    ${SRC_INCDIRS} ${SRC_DEFS} ${DEFS} $<)
	@cp ${@:%$(suffix $@)=%.d} ${@:%$(suffix $@)=%.P}; \
	 sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	     -e '/^$$/ d' -e 's/$$/ :/' < ${@:%$(suffix $@)=%.d} \
	     >> ${@:%$(suffix $@)=%.P}; \
	 rm -f ${@:%$(suffix $@)=%.d}
endef

# INCLUDE_SUBMAKEFILE - Parameterized "function" that includes a new
#   "submakefile" fragment into the overall Makefile. It also recursively
#   includes all submakefiles of the specified submakefile fragment.
#
#   USE WITH EVAL
#
define INCLUDE_SUBMAKEFILE
    # Initialize all variables that can be defined by a makefile fragment, then
    # include the specified makefile fragment.
    TARGET        :=
    TGT_CFLAGS    :=
    TGT_CXXFLAGS  :=
    TGT_DEFS      :=
    TGT_INCDIRS   :=
    TGT_LDFLAGS   :=
    TGT_LDLIBS    :=
    TGT_LINKER    :=
    TGT_POSTCLEAN :=
    TGT_POSTMAKE  :=
    TGT_PREREQS   :=

    SOURCES       :=
    SRC_CFLAGS    :=
    SRC_CXXFLAGS  :=
    SRC_DEFS      :=
    SRC_INCDIRS   :=

    SUBMAKEFILES  :=

    # A directory stack is maintained so that the correct paths are used as we
    # recursively include all submakefiles. Get the makefile's directory and
    # push it onto the stack.
    DIR := $(call CANONICAL_PATH,$(dir ${1}))
    DIR_STACK := $$(call PUSH,$${DIR_STACK},$${DIR})

    include ${1}

    # Initialize internal local variables.
    OBJS :=

    # Determine which target this makefile's variables apply to. A stack is
    # used to keep track of which target is the "current" target as we
    # recursively include other submakefiles.
    ifneq "$$(strip $${TARGET})" ""
        # This makefile defined a new target. Target variables defined by this
        # makefile apply to this new target. Initialize the target's variables.

        ifeq "$$(suffix $${TARGET})" ".a"
          TGT := $$(addprefix lib/, $$(strip $${TARGET}))
        else
          TGT := $$(addprefix bin/, $$(strip $${TARGET}))
        endif
        ALL_TGTS += $${TGT}
        $${TGT}_CFLAGS    := $${TGT_CFLAGS}
        $${TGT}_CXXFLAGS  := $${TGT_CXXFLAGS}
        $${TGT}_DEFS      := $${TGT_DEFS}
        $${TGT}_DEPS      :=
        TGT_INCDIRS       := $$(call QUALIFY_PATH,$${DIR},$${TGT_INCDIRS})
        TGT_INCDIRS       := $$(call CANONICAL_PATH,$${TGT_INCDIRS})
        $${TGT}_INCDIRS   := $${TGT_INCDIRS}
        $${TGT}_LDFLAGS   := $${TGT_LDFLAGS}
        $${TGT}_LDLIBS    := $${TGT_LDLIBS}
        $${TGT}_LINKER    := $${TGT_LINKER}
        $${TGT}_OBJS      :=
        $${TGT}_POSTCLEAN := $${TGT_POSTCLEAN}
        $${TGT}_POSTMAKE  := $${TGT_POSTMAKE}
        $${TGT}_PREREQS   := $$(addprefix $${TARGET_DIR}/lib/,$${TGT_PREREQS})
        $${TGT}_SOURCES   :=
    else
        # The values defined by this makefile apply to the the "current" target
        # as determined by which target is at the top of the stack.
        TGT := $$(strip $$(call PEEK,$${TGT_STACK}))
        $${TGT}_CFLAGS    += $${TGT_CFLAGS}
        $${TGT}_CXXFLAGS  += $${TGT_CXXFLAGS}
        $${TGT}_DEFS      += $${TGT_DEFS}
        TGT_INCDIRS       := $$(call QUALIFY_PATH,$${DIR},$${TGT_INCDIRS})
        TGT_INCDIRS       := $$(call CANONICAL_PATH,$${TGT_INCDIRS})
        $${TGT}_INCDIRS   += $${TGT_INCDIRS}
        $${TGT}_LDFLAGS   += $${TGT_LDFLAGS}
        $${TGT}_LDLIBS    += $${TGT_LDLIBS}
        $${TGT}_POSTCLEAN += $${TGT_POSTCLEAN}
        $${TGT}_POSTMAKE  += $${TGT_POSTMAKE}
        $${TGT}_PREREQS   += $${TGT_PREREQS}
    endif

    # Push the current target onto the target stack.
    TGT_STACK := $$(call PUSH,$${TGT_STACK},$${TGT})

    ifneq "$$(strip $${SOURCES})" ""
        # This makefile builds one or more objects from source. Validate the
        # specified sources against the supported source file types.
        BAD_SRCS := $$(strip $$(filter-out $${ALL_SRC_EXTS},$${SOURCES}))
        ifneq "$${BAD_SRCS}" ""
            $$(error Unsupported source file(s) found in ${1} [$${BAD_SRCS}])
        endif

        # Qualify and canonicalize paths.
        SOURCES     := $$(call QUALIFY_PATH,$${DIR},$${SOURCES})
        SOURCES     := $$(call CANONICAL_PATH,$${SOURCES})
        SRC_INCDIRS := $$(call QUALIFY_PATH,$${DIR},$${SRC_INCDIRS})
        SRC_INCDIRS := $$(call CANONICAL_PATH,$${SRC_INCDIRS})

        # Save the list of source files for this target.
        $${TGT}_SOURCES += $${SOURCES}

        # Convert the source file names to their corresponding object file
        # names.
        OBJS := $$(addprefix $${BUILD_DIR}/$$(call CANONICAL_PATH,$${TGT})/,\
                   $$(addsuffix .o,$$(basename $${SOURCES})))

        # Add the objects to the current target's list of objects, and create
        # target-specific variables for the objects based on any source
        # variables that were defined.
        $${TGT}_OBJS += $${OBJS}
        $${TGT}_DEPS += $${OBJS:%.o=%.P}
        $${OBJS}: SRC_CFLAGS   := $${$${TGT}_CFLAGS} $${SRC_CFLAGS}
        $${OBJS}: SRC_CXXFLAGS := $${$${TGT}_CXXFLAGS} $${SRC_CXXFLAGS}
        $${OBJS}: SRC_DEFS     := $$(addprefix -D,$${$${TGT}_DEFS} $${SRC_DEFS})
        $${OBJS}: SRC_INCDIRS  := $$(addprefix -I,\
                                     $${$${TGT}_INCDIRS} $${SRC_INCDIRS})
    endif

    ifneq "$$(strip $${SUBMAKEFILES})" ""
        # This makefile has submakefiles. Recursively include them.
        $$(foreach MK,$${SUBMAKEFILES},\
           $$(eval $$(call INCLUDE_SUBMAKEFILE,\
                      $$(call CANONICAL_PATH,\
                         $$(call QUALIFY_PATH,$${DIR},$${MK})))))
    endif

    # Reset the "current" target to it's previous value.
    TGT_STACK := $$(call POP,$${TGT_STACK})
    TGT := $$(call PEEK,$${TGT_STACK})

    # Reset the "current" directory to it's previous value.
    DIR_STACK := $$(call POP,$${DIR_STACK})
    DIR := $$(call PEEK,$${DIR_STACK})
endef


# MIN - Parameterized "function" that results in the minimum lexical value of
#   the two values given.
define MIN
$(firstword $(sort ${1} ${2}))
endef

# PEEK - Parameterized "function" that results in the value at the top of the
#   specified colon-delimited stack.
define PEEK
$(lastword $(subst :, ,${1}))
endef

# POP - Parameterized "function" that pops the top value off of the specified
#   colon-delimited stack, and results in the new value of the stack. Note that
#   the popped value cannot be obtained using this function; use peek for that.
define POP
${1:%:$(lastword $(subst :, ,${1}))=%}
endef

# PUSH - Parameterized "function" that pushes a value onto the specified colon-
#   delimited stack, and results in the new value of the stack.
define PUSH
${2:%=${1}:%}
endef

# QUALIFY_PATH - Given a "root" directory and one or more paths, qualifies the
#   paths using the "root" directory (i.e. appends the root directory name to
#   the paths) except for paths that are absolute.
define QUALIFY_PATH
$(addprefix ${1}/,$(filter-out /%,${2})) $(filter /%,${2})
endef

###############################################################################
#
# Start of Makefile Evaluation
#
###############################################################################

# Older versions of GNU Make lack capabilities needed by boilermake.
# With older versions, "make" may simply output "nothing to do", likely leading
# to confusion. To avoid this, check the version of GNU make up-front and
# inform the user if their version of make doesn't meet the minimum required.

MIN_MAKE_VERSION := 3.81
MIN_MAKE_VER_MSG := boilermake requires GNU Make ${MIN_MAKE_VERSION} or greater
ifeq "${MAKE_VERSION}" ""
    $(info GNU Make not detected)
    $(error ${MIN_MAKE_VER_MSG})
endif
ifneq "${MIN_MAKE_VERSION}" "$(call MIN,${MIN_MAKE_VERSION},${MAKE_VERSION})"
    $(info This is GNU Make version ${MAKE_VERSION})
    $(error ${MIN_MAKE_VER_MSG})
endif

# Define the source file extensions that we know how to handle.

C_SRC_EXTS := %.c
CXX_SRC_EXTS := %.C %.cc %.cp %.cpp %.CPP %.cxx %.c++
JAVA_EXTS    := %.jar %.tar
ALL_SRC_EXTS := ${C_SRC_EXTS} ${CXX_SRC_EXTS} ${JAVA_EXTS}

# Initialize global variables.

ALL_TGTS :=
DEFS :=
DIR_STACK :=
INCDIRS :=
TGT_STACK :=

# Discover our OS and architecture.  These were previously used to set
# BUILD_DIR and TARGET_DIR to allow multi-platform builds.  DESTDIR will do
# that for us too.

OSTYPE      := $(shell echo `uname`)
OSVERSION   := $(shell echo `uname -r`)
MACHINETYPE := $(shell echo `uname -m`)

ifeq (${MACHINETYPE}, x86_64)
  MACHINETYPE = amd64
endif

ifeq (${MACHINETYPE}, Power Macintosh)
  MACHINETYPE = ppc
endif

ifeq (${OSTYPE}, SunOS)
  MACHINETYPE = ${shell echo `uname -p`}
  ifeq (${MACHINETYPE}, sparc)
    ifeq (${shell /usr/bin/isainfo -b}, 64)
      MACHINETYPE = sparc64
    else
      MACHINETYPE = sparc32
    endif
  endif
endif

#  Set paths for building and installing.  If DESTDIR doesn't exist, use the
#  directory just above us.

ifeq "$(strip ${DESTDIR})" ""
  BUILD_DIR    := $(realpath ..)/build/obj
  TARGET_DIR   := $(realpath ..)/build
  OBSOLETE_DIR := $(realpath ..)/${OSTYPE}-${MACHINETYPE}
else
  BUILD_DIR    := $(DESTDIR)/$(MODULE)/build/obj
  TARGET_DIR   := $(DESTDIR)/$(MODULE)/build
  OBSOLETE_DIR := $(DESTDIR)/$(MODULE)/${OSTYPE}-${MACHINETYPE}
endif

#
#  Set compiler and flags based on discovered hardware
#
#  By default, debug symbols are included in all builds (even optimized).
#
#  BUILDOPTIMIZED  will disable debug symbols (leaving it just optimized)
#  BUILDDEBUG      will disable optimization  (leaving it just with debug symbols)
#  BUILDSTACKTRACE will enable stack trace on crashes, only works for Linux
#                  set to 0 on command line to disable (it's enabled by default for Linux)
#
#  BUILDPROFILE used to add -pg to LDFLAGS, and remove -D_GLIBCXX_PARALLE from CXXFLAGS and LDFLAGS,
#  and remove -fomit-frame-pointer from CXXFLAGS.  It added a bunch of complication and wasn't
#  really used.
#
#  BUILDJEMALLOC will enable jemalloc library support.
#

ifeq ($(origin CXXFLAGS), undefined)
  ifeq ($(BUILDOPTIMIZED), 1)
  else
    CXXFLAGS += -g3
  endif

  ifeq ($(BUILDDEBUG), 1)
  else
    CXXFLAGS += -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
  endif

  ifeq ($(BUILDJEMALLOC), 1)
    CXXFLAGS += -DJEMALLOC -I`jemalloc-config --includedir`
    LDFLAGS  += -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -ljemalloc `jemalloc-config --libs`
  endif

  #  Enable some warnings.
  #     gcc7:  -Wno-format-truncation  - sprintf() into possibly smaller buffer
  #            -Wno-parentheses
  CXXFLAGS += -Wall -Wextra -Wformat
  CXXFLAGS += -Wno-char-subscripts
  CXXFLAGS += -Wno-sign-compare
  CXXFLAGS += -Wno-unused-function
  CXXFLAGS += -Wno-unused-parameter
  CXXFLAGS += -Wno-unused-variable
  CXXFLAGS += -Wno-deprecated-declarations
  CXXFLAGS += -Wno-format-truncation
  CXXFLAGS += -std=c++11
else
  CXXFLAGSUSER := ${CXXFLAGS}
endif




ifeq (${OSTYPE}, Linux)
  CC        ?= gcc
  CXX       ?= g++

  CXXFLAGS  += -pthread -fopenmp -fPIC
  LDFLAGS   += -pthread -fopenmp -lm

  BUILDSTACKTRACE ?= 1
endif


#  The default MacOS compiler - even as of 10.13 High Sierra - doesn't support OpenMP.
#  Clang 6.0 installed from MacPorts supports OpenMP, but fails to compile Canu.
#  So, we require gcc7 (from MacPorts) or gcc8 (from hommebrew).
#
#  If from MacPorts:
#    port install gcc7
#    port select gcc mp-gcc7
#
#  If CC is set to 'cc', the GNU make default, we'll automagically search for other
#  versions and use those if found, preferring gcc7 over gcc8.
#
#  There' definitely a clever way to do this with 'foreach', but my Make is lacking.
#
ifeq (${OSTYPE}, Darwin)
  ifeq ($(CC), cc)
    CC7    := $(shell echo `which gcc-mp-7`)
    CXX7   := $(shell echo `which g++-mp-7`)

    ifdef CXX7
      CC  := $(CC7)
      CXX := $(CXX7)
    endif
  endif

  ifeq ($(CC), cc)
    CC8    := $(shell echo `which gcc-mp-8`)
    CXX8   := $(shell echo `which g++-mp-8`)

    ifdef CXX8
      CC  := $(CC8)
      CXX := $(CXX8)
    endif
  endif

  ifeq ($(CC), cc)
    CC9    := $(shell echo `which gcc-mp-9`)
    CXX9   := $(shell echo `which g++-mp-9`)

    ifdef CXX9
      CC  := $(CC9)
      CXX := $(CXX9)
    endif
  endif

  ifeq ($(CC), cc)
    CC8    := $(shell echo `which gcc-7`)
	  CXX8   := $(shell echo `which g++-7`)

    ifdef CXX8
      CC  := $(CC8)
      CXX := $(CXX8)
    endif
  endif

  ifeq ($(CC), cc)
    CC8    := $(shell echo `which gcc-8`)
	  CXX8   := $(shell echo `which g++-8`)

    ifdef CXX8
      CC  := $(CC8)
      CXX := $(CXX8)
    endif
  endif

  ifeq ($(CC), cc)
    CC9    := $(shell echo `which gcc-9`)
	  CXX9   := $(shell echo `which g++-9`)

    ifdef CXX9
      CC  := $(CC9)
      CXX := $(CXX9)
    endif
  endif

  ifneq ($(shell echo `$(CXX) --version 2>&1 | grep -c clang`), 0)
     CPATH := $(shell echo `which $(CXX)`)
     CLANG := $(shell echo `$(CXX) --version 2>&1 | grep clang`)
     space := 

     $(warning )
     ifeq ($(CXX), $(CPATH))
       $(warning Compiler '$(CXX)' reports version '$(CLANG)'.)
     else
       $(warning Compiler '$(CXX)' at '$(CPATH)' reports version '$(CLANG)'.)
     endif
     $(warning )
     $(warning Canu cannot be compiled with this compiler.  Please install GCC and/or)
     $(warning specify a non-Clang compiler on the command line, e.g.,)   #  Quite the evil trick to get
     $(warning $(space)    make CC=/path/to/gcc CXX=/path/to/g++);        #  this line indented!
     $(warning )
     $(error unsupported compiler)
  endif

  CXXFLAGS += -fopenmp -pthread -fPIC -m64 -Wno-format
  LDFLAGS  += -fopenmp -pthread -lm
endif


ifeq (${OSTYPE}, FreeBSD)
ifeq (${CANU_BUILD_ENV}, ports)

  # If building in the FreeBSD ports system, use the architecture as defined
  # there (technically, -p, not -m) and assume compiler and most options
  # are already defined correctly.

  MACHINETYPE=${ARCH}

  CXXFLAGS  += -pthread -fopenmp -fPIC
  LDFLAGS   += -pthread -fopenmp

else

  # Ignore the gmake default 'c++' and force g++9.
	ifeq ($(origin CXX), default)
    CC    = gcc9
    CXX   = g++9
    CCLIB = -rpath /usr/local/lib/gcc9
  endif

  #  GCC
  CXXFLAGS  += -I/usr/local/include -pthread -fopenmp -fPIC
  LDFLAGS   += -L/usr/local/lib     -pthread -fopenmp ${CCLIB} -lm -lexecinfo

  #  CLANG
  #CXXFLAGS  += -I/usr/local/include -pthread -fPIC
  #LDFLAGS   += -L/usr/local/lib     -pthread -lm -lexecinfo -lgomp

  #  Google Performance Tools malloc and heapchecker (HEAPCHECK=normal)
  #CXXFLAGS  +=
  #LDFLAGS   += -ltcmalloc

  #  Google Performance Tools cpu profiler (CPUPROFILE=/path)
  #CXXFLAGS  +=
  #LDFLAGS   += -lprofiler

  #  callgrind
  #CXXFLAGS  += -g3 -Wa,--gstabs -save-temps
endif
endif


ifneq (,$(findstring CYGWIN, ${OSTYPE}))
  CC        ?= gcc
  CXX       ?= g++

  CXXFLAGS  := -fopenmp -pthread
  LDFLAGS   := -fopenmp -pthread -lm
endif


#  Stack tracing support.  Wow, what a pain.  Only Linux is supported.  This is just documentation,
#  don't actually enable any of this stuff!
#
#  backward-cpp looks very nice, only a single header file.  But it needs libberty (-liberty) and
#  libbfd (-lbfd).  The former should be installed with gcc, and the latter is in elfutils.  On
#  three out of our three development machines, it fails for various reasons.
#
#  libunwind is pretty basic.
#
#  libbacktrace works (on Linux) and is simple enough to include in our tree.
#
#  None of these give any useful information on BSDs (which includes OS X aka macOS).
#
#
#  Backtraces with libunwind.  Not informative on FreeBSD.
#CXXFLAGS  += -DLIBUNWIND
#LDFLAGS   +=
#LDLIBS    += -lunwind -lunwind-x86_64
#
#
#  Backtraces with libbacktrace.  FreeBSD works, but trace is empty.
#BUILDSTACK = 1
#CXXFLAGS  += -DLIBBACKTRACE
#LDFLAGS   +=
#LDLIBS    +=
#
#
#  Backtraces with backward-cpp.
#
#  Stack walking:
#    BACKWARD_HAS_UNWIND    - used by gcc/clang for exception handling
#    BACKWARD_HAS_BACKTRACE - part of glib, not as accurate, more portable
#
#  Stack interpretation:
#    BACKWARE_HAS_DW               - most information, libdw, (elfutils or libdwarf)
#    BACKWARD_HAS_BFD              - some information, libbfd
#    BACKWARD_HAS_BACKTRACE_SYMBOL - minimal information (file and function), portable
#
#  helix   fails with: cannot find -liberty
#  gryphon fails with: cannot find -lbfd
#  freebsd can't install a working elfutils, needed for libdw"
#    In file included from AS_UTL/AS_UTL_stackTrace.C:183:0:
#    AS_UTL/backward.hpp:241:30: fatal error: elfutils/libdw.h: No such file or directory
#     #  include <elfutils/libdw.h>
#
#CXXFLAGS  += -DBACKWARDCPP -DBACKWARD_HAS_BFD
#LDFLAGS   +=
#LDLIBS    += -lbfd -liberty -ldl -lz
#
#  Needs libdw, elfutils
#CXXFLAGS  += -DBACKWARDCPP -DBACKWARD_HAS_DW
#LDFLAGS   +=
#LDLIBS    += -ldl -lz
#
#  Generates nothing useful, no function names, just binary names
#CXXFLAGS  += -DBACKWARDCPP
#LDFLAGS   +=
#LDLIBS    += -ldl -lz
#
#
#  No backtrace support.
#CXXFLAGS   += -DNOBACKTRACE

#  But, if we we have an old GCC, stack tracing support isn't there.  The
#  second test is because gcc7 and gcc8 report the major version for
#  -dumpversion.

GXX_45 := $(shell expr `${CXX} -dumpversion     | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40500)
GXX_VV := $(shell ${CXX} -dumpversion)

ifneq (,$(filter ${GXX_VV}, 7 8 9))
GXX_45 := $(shell expr `${CXX} -dumpfullversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40500)
GXX_VV := $(shell ${CXX} -dumpfullversion)
endif

ifeq (${BUILDSTACKTRACE}, 1)
ifeq (${GXX_45}, 0)
$(info WARNING:)
$(info WARNING: ${CXX} ${GXX_VV} detected, disabling stack trace support.  Please upgrade to GCC 4.7 or higher.)
$(info WARNING:)
BUILDSTACKTRACE = 0
endif
endif

ifeq (${BUILDSTACKTRACE}, 1)
CXXFLAGS  += -DLIBBACKTRACE
else
CXXFLAGS  += -DNOBACKTRACE
endif


# Include the main user-supplied submakefile. This also recursively includes
# all other user-supplied submakefiles.
$(eval $(call INCLUDE_SUBMAKEFILE,main.mk))

# Perform post-processing on global variables as needed.
DEFS := $(addprefix -D,${DEFS})
INCDIRS := $(addprefix -I,$(call CANONICAL_PATH,${INCDIRS}))

# Define the "all" target (which simply builds all user-defined targets) as the
# default goal.
.PHONY: all
all: $(addprefix ${TARGET_DIR}/,${ALL_TGTS}) \
     ${TARGET_DIR}/bin/canu \
     ${TARGET_DIR}/bin/canu-time \
     ${TARGET_DIR}/bin/draw-tig \
     ${TARGET_DIR}/bin/canu.defaults \
     ${TARGET_DIR}/share/java/classes/mhap-2.1.3.jar \
     ${TARGET_DIR}/lib/site_perl/canu/Consensus.pm \
     ${TARGET_DIR}/lib/site_perl/canu/CorrectReads.pm \
     ${TARGET_DIR}/lib/site_perl/canu/HaplotypeReads.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Configure.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Defaults.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Execution.pm \
     ${TARGET_DIR}/lib/site_perl/canu/SequenceStore.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Grid.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Grid_Cloud.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Grid_DNANexus.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Grid_LSF.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Grid_PBSTorque.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Grid_SGE.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Grid_Slurm.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Meryl.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Output.pm \
     ${TARGET_DIR}/lib/site_perl/canu/OverlapBasedTrimming.pm \
     ${TARGET_DIR}/lib/site_perl/canu/OverlapErrorAdjustment.pm \
     ${TARGET_DIR}/lib/site_perl/canu/OverlapInCore.pm \
     ${TARGET_DIR}/lib/site_perl/canu/OverlapMhap.pm \
     ${TARGET_DIR}/lib/site_perl/canu/OverlapMMap.pm \
     ${TARGET_DIR}/lib/site_perl/canu/OverlapStore.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Report.pm \
     ${TARGET_DIR}/lib/site_perl/canu/Unitig.pm
	@if [ -e ${OBSOLETE_DIR}/bin/canu ] ; then echo "Moving obsolete installation from ${OBSOLETE_DIR} to ${OBSOLETE_DIR}-OBSOLETE."; mv ${OBSOLETE_DIR} ${OBSOLETE_DIR}-OBSOLETE ; fi
	@echo ""
	@echo "Success!"
	@echo "canu installed in ${TARGET_DIR}/bin/canu"
	@echo ""

# Add a new target rule for each user-defined target.
$(foreach TGT,${ALL_TGTS},\
  $(eval $(call ADD_TARGET_RULE,${TGT})))

# Add pattern rule(s) for creating compiled object code from C source.
$(foreach TGT,${ALL_TGTS},\
  $(foreach EXT,${C_SRC_EXTS},\
    $(eval $(call ADD_OBJECT_RULE,${BUILD_DIR}/$(call CANONICAL_PATH,${TGT}),\
             ${EXT},$${COMPILE_C_CMDS}))))

# Add pattern rule(s) for creating compiled object code from C++ source.
$(foreach TGT,${ALL_TGTS},\
  $(foreach EXT,${CXX_SRC_EXTS},\
    $(eval $(call ADD_OBJECT_RULE,${BUILD_DIR}/$(call CANONICAL_PATH,${TGT}),\
             ${EXT},$${COMPILE_CXX_CMDS}))))

# Add "clean" rules to remove all build-generated files.
.PHONY: clean
$(foreach TGT,${ALL_TGTS},\
  $(eval $(call ADD_CLEAN_RULE,${TGT})))

# Include generated rules that define additional (header) dependencies.
$(foreach TGT,${ALL_TGTS},\
  $(eval -include ${${TGT}_DEPS}))


${TARGET_DIR}/bin/canu: pipelines/canu.pl
	cp -pf pipelines/canu.pl ${TARGET_DIR}/bin/canu
	@chmod +x ${TARGET_DIR}/bin/canu

${TARGET_DIR}/bin/canu-time: pipelines/canu-time.pl
	cp -pf pipelines/canu-time.pl ${TARGET_DIR}/bin/canu-time
	@chmod +x ${TARGET_DIR}/bin/canu-time

${TARGET_DIR}/bin/draw-tig: pipelines/draw-tig.pl
	cp -pf pipelines/draw-tig.pl ${TARGET_DIR}/bin/draw-tig
	@chmod +x ${TARGET_DIR}/bin/draw-tig


${TARGET_DIR}/bin/canu.defaults:
	@echo > ${TARGET_DIR}/bin/canu.defaults  "# Add site specific options (for setting up Grid or limiting memory/threads) here."


${TARGET_DIR}/share/java/classes/mhap-2.1.3.jar: mhap/mhap-2.1.3.jar
	cp -pf mhap/mhap-2.1.3.jar ${TARGET_DIR}/share/java/classes/mhap-2.1.3.jar


${TARGET_DIR}/lib/site_perl/canu/Consensus.pm: pipelines/canu/Consensus.pm
	cp -pf pipelines/canu/Consensus.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/CorrectReads.pm: pipelines/canu/CorrectReads.pm
	cp -pf pipelines/canu/CorrectReads.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/HaplotypeReads.pm: pipelines/canu/HaplotypeReads.pm
	cp -pf pipelines/canu/HaplotypeReads.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Configure.pm: pipelines/canu/Configure.pm
	cp -pf pipelines/canu/Configure.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Defaults.pm: pipelines/canu/Defaults.pm
	cp -pf pipelines/canu/Defaults.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/ErrorEstimate.pm: pipelines/canu/ErrorEstimate.pm
	cp -pf pipelines/canu/ErrorEstimate.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Execution.pm: pipelines/canu/Execution.pm
	cp -pf pipelines/canu/Execution.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/SequenceStore.pm: pipelines/canu/SequenceStore.pm
	cp -pf pipelines/canu/SequenceStore.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Grid.pm: pipelines/canu/Grid.pm
	cp -pf pipelines/canu/Grid.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Grid_Cloud.pm: pipelines/canu/Grid_Cloud.pm
	cp -pf pipelines/canu/Grid_Cloud.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Grid_DNANexus.pm: pipelines/canu/Grid_DNANexus.pm
	cp -pf pipelines/canu/Grid_DNANexus.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Grid_LSF.pm: pipelines/canu/Grid_LSF.pm
	cp -pf pipelines/canu/Grid_LSF.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Grid_PBSTorque.pm: pipelines/canu/Grid_PBSTorque.pm
	cp -pf pipelines/canu/Grid_PBSTorque.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Grid_SGE.pm: pipelines/canu/Grid_SGE.pm
	cp -pf pipelines/canu/Grid_SGE.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Grid_Slurm.pm: pipelines/canu/Grid_Slurm.pm
	cp -pf pipelines/canu/Grid_Slurm.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Meryl.pm: pipelines/canu/Meryl.pm
	cp -pf pipelines/canu/Meryl.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Output.pm: pipelines/canu/Output.pm
	cp -pf pipelines/canu/Output.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/OverlapBasedTrimming.pm: pipelines/canu/OverlapBasedTrimming.pm
	cp -pf pipelines/canu/OverlapBasedTrimming.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/OverlapErrorAdjustment.pm: pipelines/canu/OverlapErrorAdjustment.pm
	cp -pf pipelines/canu/OverlapErrorAdjustment.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/OverlapInCore.pm: pipelines/canu/OverlapInCore.pm
	cp -pf pipelines/canu/OverlapInCore.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/OverlapMhap.pm: pipelines/canu/OverlapMhap.pm
	cp -pf pipelines/canu/OverlapMhap.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/OverlapMMap.pm: pipelines/canu/OverlapMMap.pm
	cp -pf pipelines/canu/OverlapMMap.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/OverlapStore.pm: pipelines/canu/OverlapStore.pm
	cp -pf pipelines/canu/OverlapStore.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Report.pm: pipelines/canu/Report.pm
	cp -pf pipelines/canu/Report.pm ${TARGET_DIR}/lib/site_perl/canu/

${TARGET_DIR}/lib/site_perl/canu/Unitig.pm: pipelines/canu/Unitig.pm
	cp -pf pipelines/canu/Unitig.pm ${TARGET_DIR}/lib/site_perl/canu/

#  Makefile processed.  Regenerate the version number file, make some
#  directories, and report that we're starting the build.

$(eval $(shell ../scripts/version_update.pl $(MODULE) utility/src/utility/version.H))

$(shell mkdir -p ${TARGET_DIR}/lib/site_perl/canu)
$(shell mkdir -p ${TARGET_DIR}/share/java/classes)
$(shell mkdir -p ${TARGET_DIR}/bin)

$(info For '${OSTYPE}' '${OSVERSION}' as '${MACHINETYPE}' into '${TARGET_DIR}/{bin,obj}'.)
$(info Using '$(shell which ${CXX})' version '${GXX_VV}'.)
ifneq ($(origin CXXFLAGSUSER), undefined)
$(info Using user-supplied CXXFLAGS '${CXXFLAGSUSER}'.)
endif
$(info )
