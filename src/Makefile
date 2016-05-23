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
${1}/%.o: ${2}
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
    ifeq "$$(suffix ${1})" ".jar"
        # Add a target for copying files.  Actually, just untars the files, kind
        # of ugly, but the only way I could discover to get the lib/ files installed.
        $${TARGET_DIR}/${1}: $${${1}}
	    @mkdir -p $$(dir $$@)
	    $$(strip tar -C $$(dir $$@) -xf $${${1}_SOURCES})
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
	    $$(strip $${${1}_LINKER} -o $$@ $${LDFLAGS} $${${1}_LDFLAGS} \
	        $${${1}_OBJS} $${LDLIBS} $${${1}_LDLIBS})
	    $${${1}_POSTMAKE}
    endif
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
        TGT := $$(strip $${TARGET})
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
        $${TGT}_PREREQS   := $$(addprefix $${TARGET_DIR}/,$${TGT_PREREQS})
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

# Discover our OS and architecture.  These are used to set the BUILD_DIR and TARGET_DIR to
# something more useful than 'build' and '.'.

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

# Set compiler and flags based on discovered hardware

ifeq (${OSTYPE}, Linux)
  CC        ?= gcc
  CXX       ?= g++
  CXXFLAGS  := -pthread -fPIC -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare # -Wformat -Wall -Wextra
  LDFLAGS   := -pthread -lm

  CXXFLAGS  += -D_GLIBCXX_PARALLEL -fopenmp
  LDFLAGS   +=                     -fopenmp

  CXXFLAGS  += -Wall -Wextra -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -Wformat

  ifeq ($(BUILDPROFILE), 1)
    CXXFLAGS  +=
    LDFLAGS   += -pg
  else
    CXXFLAGS  += -D_GLIBCXX_PARALLEL
    LDFLAGS   += -D_GLIBCXX_PARALLEL
  endif

  ifeq ($(BUILDDEBUG), 1)
    CXXFLAGS  += -g
  else
    ifeq ($(BUILDPROFILE), 1)
      CXXFLAGS   += -O4 -funroll-loops -fexpensive-optimizations -finline-functions
    else
      CXXFLAGS   += -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
    endif
  endif
endif


ifeq (${OSTYPE}, Darwin)
  CLANG = $(shell echo `${CXX} --version 2>&1 | grep -c clang`)

  CC       ?= gcc
  CXX      ?= g++
  ifeq ($(CLANG), 0)
     CXXFLAGS := -pthread -fopenmp -fPIC -m64 -Wall -Wextra -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -Wformat
     LDFLAGS  := -pthread -fopenmp -lm
  else
     CXXFLAGS := -DBROKEN_CLANG_OpenMP -pthread -fPIC -m64 -Wall -Wextra -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -Wformat
     LDFLAGS  := -pthread -lm
  endif

  ifeq ($(BUILDPROFILE), 1)
    CXXFLAGS  +=
    LDFLAGS   += -pg
  else
    CXXFLAGS  += -D_GLIBCXX_PARALLEL
    LDFLAGS   += -D_GLIBCXX_PARALLEL
  endif

  CXXFLAGS  +=
  LDFLAGS   +=

  ifeq ($(BUILDDEBUG), 1)
    CXXFLAGS  += -g
# -fno-omit-frame-pointer
  else
    ifeq ($(BUILDPROFILE), 1)
      CXXFLAGS   += -g3 -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fno-omit-frame-pointer
    else
      CXXFLAGS   +=     -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
    endif
  endif
endif


ifeq (${OSTYPE}, FreeBSD)
ifeq (${MACHINETYPE}, amd64)
  CC       ?= gcc48
  CXX      ?= g++48

  CXXFLAGS := -pthread -fPIC -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -Wformat -Wextra
  LDFLAGS  := -pthread -lm

  CXXFLAGS  += -fopenmp
  LDFLAGS   += -fopenmp -rpath /usr/local/lib/gcc48

  #  Google Performance Tools malloc and heapchecker (HEAPCHECK=normal)
  #CXXFLAGS  +=
  #LDFLAGS   += -ltcmalloc

  #  Google Performance Tools cpu profiler (CPUPROFILE=/path)
  #CXXFLAGS  +=
  #LDFLAGS   += -lprofiler

  ifeq ($(BUILDPROFILE), 1)
    CXXFLAGS  +=
    LDFLAGS   += -pg
  else
    CXXFLAGS  += -D_GLIBCXX_PARALLEL
    LDFLAGS   += -D_GLIBCXX_PARALLEL
  endif

  CXXFLAGS  +=             -I/usr/local/include
  LDFLAGS   +=             -L/usr/local/lib

  #  callgrind
  #CXXFLAGS  += -g3 -Wa,--gstabs -save-temps

  #CXXFLAGS  += -DLIBUNWIND -I/usr/local/include
  #LDFLAGS   +=             -L/usr/local/lib -lunwind -lunwind-x86_64

  ifeq ($(BUILDDEBUG), 1)
    CXXFLAGS  += -g
  else
    ifeq ($(BUILDPROFILE), 1)
      CXXFLAGS   += -g3 -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fno-omit-frame-pointer
    else
      CXXFLAGS   +=     -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
    endif
  endif
endif
endif


ifeq (${OSTYPE}, FreeBSD)
ifeq (${MACHINETYPE}, arm)
  CC       ?= gcc48
  CXX      ?= g++48

  CXXFLAGS  := -pthread -fPIC -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -Wformat -Wextra
  LDFLAGS   := -pthread -lm

  CXXFLAGS  += -D_GLIBCXX_PARALLEL -fopenmp
  LDFLAGS   += -D_GLIBCXX_PARALLEL -fopenmp -rpath /usr/local/lib/gcc48

  CXXFLAGS  += -I/usr/local/include -O3 -funroll-loops -fomit-frame-pointer
  LDFLAGS   += -L/usr/local/lib
endif
endif

ifneq (,$(findstring CYGWIN, ${OSTYPE}))
  CC        ?= gcc
  CXX       ?= g++
  CXXFLAGS  := -pthread -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare
  LDFLAGS   := -pthread -lm

  CXXFLAGS  += -fopenmp
  LDFLAGS   += -fopenmp

  CXXFLAGS  += -Wall -Wextra -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -Wformat

  ifeq ($(BUILDPROFILE), 1)
    CXXFLAGS  +=
    LDFLAGS   += -pg
  else
    CXXFLAGS  += 
    LDFLAGS   += 
  endif

  ifeq ($(BUILDDEBUG), 1)
    CXXFLAGS  += -g
  else
    ifeq ($(BUILDPROFILE), 1)
      CXXFLAGS   += -O4 -funroll-loops -fexpensive-optimizations -finline-functions
    else
      CXXFLAGS   += -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
    endif
  endif
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
all: UPDATE_VERSION MAKE_DIRS \
     $(addprefix ${TARGET_DIR}/,${ALL_TGTS}) \
     ${TARGET_DIR}/canu \
     ${TARGET_DIR}/canu.defaults \
     ${TARGET_DIR}/lib/canu/lib \
     ${TARGET_DIR}/lib/canu/Consensus.pm \
     ${TARGET_DIR}/lib/canu/CorrectReads.pm \
     ${TARGET_DIR}/lib/canu/Configure.pm \
     ${TARGET_DIR}/lib/canu/Defaults.pm \
     ${TARGET_DIR}/lib/canu/ErrorEstimate.pm \
     ${TARGET_DIR}/lib/canu/Execution.pm \
     ${TARGET_DIR}/lib/canu/Gatekeeper.pm \
     ${TARGET_DIR}/lib/canu/Grid.pm \
     ${TARGET_DIR}/lib/canu/Grid_LSF.pm \
     ${TARGET_DIR}/lib/canu/Grid_PBSTorque.pm \
     ${TARGET_DIR}/lib/canu/Grid_SGE.pm \
     ${TARGET_DIR}/lib/canu/Grid_Slurm.pm \
     ${TARGET_DIR}/lib/canu/HTML.pm \
     ${TARGET_DIR}/lib/canu/Meryl.pm \
     ${TARGET_DIR}/lib/canu/Output.pm \
     ${TARGET_DIR}/lib/canu/OverlapBasedTrimming.pm \
     ${TARGET_DIR}/lib/canu/OverlapErrorAdjustment.pm \
     ${TARGET_DIR}/lib/canu/OverlapInCore.pm \
     ${TARGET_DIR}/lib/canu/OverlapMhap.pm \
     ${TARGET_DIR}/lib/canu/OverlapMMap.pm \
     ${TARGET_DIR}/lib/canu/OverlapStore.pm \
     ${TARGET_DIR}/lib/canu/Unitig.pm
	@echo ""
	@echo "Success!"
	@echo "canu installed in ${TARGET_DIR}/canu"

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

#  A fake target, to regenerate the canu_version.H file on every build.
.PHONY: UPDATE_VERSION
UPDATE_VERSION:
	@./canu_version_update.pl

AS_global.C: UPDATE_VERSION

#  A fake target, to make the directory for the canu perl modules.
.PHONY: MAKE_DIRS
MAKE_DIRS:
	@if [ ! -e ${TARGET_DIR}/lib/canu ] ; then mkdir -p ${TARGET_DIR}/lib/canu ; fi

${TARGET_DIR}/canu: pipelines/canu.pl
	cp -pf pipelines/canu.pl ${TARGET_DIR}/canu
	chmod +x ${TARGET_DIR}/canu

${TARGET_DIR}/canu.defaults:
	echo > ${TARGET_DIR}/canu.defaults  "# Add site specific options (for setting up Grid or limiting memory/threads) here."
	chmod -x ${TARGET_DIR}/canu.defaults

${TARGET_DIR}/lib/canu/lib: pipelines/install-perl-libraries.sh
	cd pipelines && sh install-perl-libraries.sh ${TARGET_DIR}/lib/canu

${TARGET_DIR}/lib/canu/Consensus.pm: pipelines/canu/Consensus.pm
	cp -pf pipelines/canu/Consensus.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/CorrectReads.pm: pipelines/canu/CorrectReads.pm
	cp -pf pipelines/canu/CorrectReads.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Configure.pm: pipelines/canu/Configure.pm
	cp -pf pipelines/canu/Configure.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Defaults.pm: pipelines/canu/Defaults.pm
	cp -pf pipelines/canu/Defaults.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/ErrorEstimate.pm: pipelines/canu/ErrorEstimate.pm
	cp -pf pipelines/canu/ErrorEstimate.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Execution.pm: pipelines/canu/Execution.pm
	cp -pf pipelines/canu/Execution.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Gatekeeper.pm: pipelines/canu/Gatekeeper.pm
	cp -pf pipelines/canu/Gatekeeper.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Grid.pm: pipelines/canu/Grid.pm
	cp -pf pipelines/canu/Grid.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Grid_LSF.pm: pipelines/canu/Grid_LSF.pm
	cp -pf pipelines/canu/Grid_LSF.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Grid_PBSTorque.pm: pipelines/canu/Grid_PBSTorque.pm
	cp -pf pipelines/canu/Grid_PBSTorque.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Grid_SGE.pm: pipelines/canu/Grid_SGE.pm
	cp -pf pipelines/canu/Grid_SGE.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Grid_Slurm.pm: pipelines/canu/Grid_Slurm.pm
	cp -pf pipelines/canu/Grid_Slurm.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/HTML.pm: pipelines/canu/HTML.pm
	cp -pf pipelines/canu/HTML.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Meryl.pm: pipelines/canu/Meryl.pm
	cp -pf pipelines/canu/Meryl.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Output.pm: pipelines/canu/Output.pm
	cp -pf pipelines/canu/Output.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/OverlapBasedTrimming.pm: pipelines/canu/OverlapBasedTrimming.pm
	cp -pf pipelines/canu/OverlapBasedTrimming.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/OverlapErrorAdjustment.pm: pipelines/canu/OverlapErrorAdjustment.pm
	cp -pf pipelines/canu/OverlapErrorAdjustment.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/OverlapInCore.pm: pipelines/canu/OverlapInCore.pm
	cp -pf pipelines/canu/OverlapInCore.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/OverlapMhap.pm: pipelines/canu/OverlapMhap.pm
	cp -pf pipelines/canu/OverlapMhap.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/OverlapMMap.pm: pipelines/canu/OverlapMMap.pm
	cp -pf pipelines/canu/OverlapMMap.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/OverlapStore.pm: pipelines/canu/OverlapStore.pm
	cp -pf pipelines/canu/OverlapStore.pm ${TARGET_DIR}/lib/canu/

${TARGET_DIR}/lib/canu/Unitig.pm: pipelines/canu/Unitig.pm
	cp -pf pipelines/canu/Unitig.pm ${TARGET_DIR}/lib/canu/

#  Makefile processed.  Report that we're starting the build.

${info Building for '${OSTYPE}' '${OSVERSION}' as '${MACHINETYPE}' into '${DESTDIR}${PREFIX}/$(OSTYPE)-$(MACHINETYPE)/{bin,obj}'}
#${info Using LD_RUN_PATH '${LD_RUN_PATH}'}
