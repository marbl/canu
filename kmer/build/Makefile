# -*- makefile -*-

### NOTE: make-3.80 has a bug in eval which causes it not to return
###   an empty target, but instead it sometimes returns very weird
###   stuff that makes make segv, I think.
### You need to get it patched at:
### http://savannah.gnu.org/bugs/?func=detailbug&bug_id=1517&group_id=71

default:   all

# kind of useful for debugging
# this target just prints out the directories traversed,
# the non-recursive includes found and a bunch of targets
debug:
	@echo DIRS ${//}
	@echo RECS ${RECURSIVE}
	@echo ALL_C_EXES ${ALL_C_EXES}
	@echo ALL_C_LIBS ${ALL_C_LIBS}
	@echo ALL_CXX_EXES ${ALL_CXX_EXES}
	@echo ALL_CXX_LIBS ${ALL_CXX_LIBS}
	@echo ALL_CXX_SHLIBS ${ALL_CXX_SHLIBS}
	@echo INSTALL ${INSTALL/}

## loading Make.compilers
# Compiler information and options.
# Always get it from the same directory as this makefile.
# Note: I should really change the name of Make.compilers to
#       Makeconf or Make.config or something like that.
ifndef MAKEFILE/
  MAKEFILE/ :=$(dir $(firstword $(MAKEFILE_LIST)))
endif

ifdef MAKE_COMPILERS
 include ${MAKE_COMPILERS}
else
 include ${MAKEFILE/}Make.compilers
endif

##### non-recursive make magic
# all directories which have been included
//           :=
# recursive directories to be filtered out of //
# and handled differently
RECURSIVE    :=
# current path (empty)
/            :=

## Include -hack
# Include is wrapped in something which will push and pop /
# properly while adding newly discovered directories to //
# and keeping track of who is who's children (.SUBS).
#
# Each directory so included has its Make.include file included.
# Those Make.include files can use $/ prepended to local names
# to prevent name pollution, and define their own subtargets.
#
define Include
 $(foreach x,$(strip ${1}),$(call Include_File,$x))
endef

define Include_File
  ifeq ($(filter ${1}.,${WITHOUT_}),)
    ifeq ($(wildcard ${1}Make.include),${1}Make.include)
      $/.SUBS +=${1}.
      // +=${1}.
      ${1}.SUBS :=
      /  :=${1}
      include ${1}Make.include
      /  :=$/
    endif
  endif

endef

##### System specific includes

## WITHOUT
# If the user specifies a WITHOUT, then those paths are not
# followed.
WITHOUT_:=$(patsubst %,%.,$(strip ${WITHOUT}))

## First Make.include inclusion
# invoke the toplevel include file.
# We use 'Include_File' instead of 'Include' since $/ is empty
$(eval $(call Include_File,$/))

#### Targets which have been declared RECURSIVE are removed
# from // and processed separately.
//            :=$(filter-out ${RECURSIVE},${//})

# RECURSIVE now holds the paths which are being done legacy style
# and // holds the paths which are going to be part of the system.

######  generic pattern rules for subgoals
# rules for C/C++
%.o: %.c
	${CC} ${CFLAGS} ${CFLAGS_COMPILE} -o $@ -c $<

%.o: %.cc
	${CXX} ${CXXFLAGS} ${CXXFLAGS_COMPILE} -o $@ -c $<

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${CXXFLAGS_COMPILE} -o $@ -c $<

%.o: %.C
	${CXX} ${CXXFLAGS} ${CXXFLAGS_COMPILE} -o $@ -c $<

### Building depends goals
# using '=' instead of ':=' for that dynamic scoping
__DEPGOALS__=

# really, only C/C++ depends goals are built here right now.
# C_SRCS and CXX_SRCS are collected together and turned into
# associated *.d dependency files.

# WISHLIST: does not propagate failure to the parent make
# for some reason.  It really should.  I think the pipe to
# sed masks the exit code.
ALL_C_SRCS       :=$(strip $(foreach x,${//},${${x:.=.C_SRCS}}))
ALL_C_DEPS       :=${ALL_C_SRCS:=.d}
${ALL_C_DEPS}:%.d:%
	@ echo "making $@"
	 dir=`echo $< | sed -e's~[^/]*$$~~'`; \
	  ${CCDEP} ${CDEPFLAGS} ${CFLAGS} $< | \
	  sed -e"/:/s!^!$${dir}!" > $@
__DEPGOALS__+=$${${1:.=.C_SRCS}:=.d}

ALL_CXX_SRCS     :=$(strip $(foreach x,${//},${${x:.=.CXX_SRCS}}))
ALL_CXX_DEPS     :=${ALL_CXX_SRCS:=.d}
${ALL_CXX_DEPS}:%.d:%
	@ echo "making $@"
	 dir=`echo $< | sed -e's~[^/]*$$~~'`; \
	  ${CXXDEP} ${CXXDEPFLAGS} ${CXXFLAGS} $< | \
	  sed -e"/:/s!^!$${dir}!" > $@
__DEPGOALS__+=$${${1:.=.CXX_SRCS}:=.d}

## include the .d files
# This contional gives us a way out if something goes horribly
# wrong...
ifeq ($(filter %-clean,${MAKECMDGOALS}),)
  include ${ALL_C_DEPS} ${ALL_CXX_DEPS}
endif

### Building subgoals
#   The goals are C/C++ programs, libs, and anything else we add
#   to the mix below.  It would not be unreasonable to allow some
#   sort of extension system whereby rules get included too.

# subgoals are done with the deferred '=' not the ':='.
# This is because we want dynamic scoping.
__SUBGOALS__=

## EXE targets
ALL_C_EXES   :=$(strip $(foreach x,${//},${${x:.=.C_EXES}}))
${ALL_C_EXES}:
	${CC} ${CLDFLAGS} -o $@ $^ ${CLIBS}
__SUBGOALS__+=$${${1:.=.C_EXES}}

ALL_CXX_EXES :=$(strip $(foreach x,${//},${${x:.=.CXX_EXES}}))
${ALL_CXX_EXES}:
	${CXX} ${CXXLDFLAGS} -o $@ $^ ${CXXLIBS}
__SUBGOALS__+=$${${1:.=.CXX_EXES}}

## LIB targets
ALL_C_LIBS   :=$(strip $(foreach x,${//},${${x:.=.C_LIBS}}))
${ALL_C_LIBS}:
	rm -f $@
	${AR} ${ARFLAGS} $@ $^
__SUBGOALS__+=$${${1:.=.C_LIBS}}

ALL_CXX_LIBS     :=$(strip $(foreach x,${//},${${x:.=.CXX_LIBS}}))
${ALL_CXX_LIBS}:
	rm -f $@
	${AR} ${ARFLAGS} $@ $^
__SUBGOALS__+=$${${1:.=.CXX_LIBS}}

## Shared targets
ALL_CXX_SHLIBS   :=$(strip $(foreach x,${//},${${x:.=.CXX_SHLIBS}}))
${ALL_CXX_SHLIBS}:
	${CXX} ${CXXFLAGS} ${SHLIB_FLAGS} -o $@ $^ ${CXXLIBS}
__SUBGOALS__+=$${${1:.=.CXX_SHLIBS}}

# any other targets and target rules we dump in here
-include ${MAKEFILE/}Make.rules

####### recursive rules, applied to all subdirectory targets
# We define the basic form of the all, clean, ... rules on a
# per-path basis (the $/.all, $..clean, ... targets).  This allows
# selective targeting.


## rules for each subtarget
#  current subtargets are :
#        all, ls, clean, real-clean, depends-clean, install
#  with all, clean, real-clean, install being required targets
#  of recursive makes.

clean:         ${RECURSIVE:.=.clean}      $/.clean
define .RULE-clean
${1:.=.clean}: $${${1:.=.SUBS}:.=.clean}
	rm -f $${${1:.=.CLEAN}} ${__SUBGOALS__}
	(cd $1 && rm -rf ${C_TMP} ${CXX_TMP})

endef
$(eval $(foreach x,${//},$(call .RULE-clean,$x)))

depends-clean:                            $/.depends-clean
${RECURSIVE:.=.depends-clean}:
define .RULE-depends-clean
${1:.=.depends-clean}: $${${1:.=.SUBS}:.=.depends-clean}
	rm -f ${1:.=Make.depends} ${__DEPGOALS__}

endef
$(eval $(foreach x,${//},$(call .RULE-depends-clean,$x)))

real-clean:    ${RECURSIVE:.=.real-clean} $/.real-clean
define .RULE-real-clean
${1:.=.real-clean}: $${${1:.=.SUBS}:.=.real-clean}
	rm -f $${${1:.=.CLEAN}} ${__SUBGOALS__}
	(cd $1 && rm -rf ${C_TMP} ${CXX_TMP})
	rm -f ${1:.=Make.depends} ${__DEPGOALS__}
	rm -f $${${1:.=.REAL-CLEAN}}

endef
$(eval $(foreach x,${//},$(call .RULE-real-clean,$x)))

all:           ${RECURSIVE:.=.all}        $/.all
define .RULE-all
${1:.=.all}: $${${1:.=.SUBS}:.=.all} ${__SUBGOALS__}

endef
$(eval $(foreach x,${//},$(call .RULE-all,$x)))

## this is kind of a horrible hack right now.  Maybe there is a
## better way to do this.
install:       ${RECURSIVE:.=.install}    $/.install
define .RULE-install
${1:.=.install}: $${${1:.=.SUBS}:.=.install} ${__SUBGOALS__} \
    $${1:.=.C_INCS} \
    $${1:.=.CXX_INCS} \
    $${1:.=.SHARES}
	@echo "========| $$@ |========"
# all exes go to bin
	files='$$(strip $${${1:.=.C_EXES}} $${${1:.=.CXX_EXES}})'; \
	if [ -n "$$$${files}" ] ; then \
	  mkdir -p ${INSTALL/}bin ; \
	  cp -p $$$${files} ${INSTALL/}bin ; \
        fi
# all libs go to lib/${.LIB/}
	files='$$(strip $${${1:.=.C_LIBS}} $${${1:.=.CXX_LIBS}} $${${1:.=.CXX_SHLIBS}})'; \
	if [ -n "$$$${files}" ] ; then \
	  mkdir -p ${INSTALL/}lib/$${${1:.=.LIB/}} ; \
	  cp -p $$$${files} ${INSTALL/}lib/$${${1:.=.LIB/}} ; \
        fi
# all includes go to include/${.INCLUDE/}
	files='$$(strip $${${1:.=.C_INCS}} $${${1:.=.CXX_INCS}})'; \
	if [ -n "$$$${files}" ] ; then \
	  mkdir -p ${INSTALL/}include/$${${1:.=.INCLUDE/}} ; \
	  cp -p $$$${files} ${INSTALL/}include/$${${1:.=.INCLUDE/}} ; \
        fi
# all share files go to share/${.SHARE/}
	files='$$(strip $${${1:.=.SHARES}})'; \
	if [ -n "$$$${files}" ] ; then \
	  mkdir -p ${INSTALL/}share/$${${1:.=.SHARE/}} ; \
	  cp -p $$$${files} ${INSTALL/}share/$${${1:.=.SHARE/}} ; \
        fi

endef
$(eval $(foreach x,${//},$(call .RULE-install,$x)))

# another sort of debugging thing.  Let's me check on which targets
# actually didn't get made, or got partially made
LSOPTS:=-l
ls:                                       $/.ls
${RECURSIVE:.=.ls}:
define .RULE-ls
${1:.=.ls}: $${${1:.=.SUBS}:.=.ls}
	@ files='$$(strip ${__SUBGOALS__})'; \
	  if [ -n "$$$${files}" ] ; then \
	    ls ${LSOPTS} $$$${files} ; \
	  fi ; exit 0;

endef
$(eval $(foreach x,${//},$(call .RULE-ls,$x)))


