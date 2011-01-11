# -*- makefile -*-

### NOTE: make-3.80 has a bug in eval which causes it not to return
###   an empty target, but instead it sometimes returns very weird
###   stuff that makes make segv, I think.
### You need to get it patched at:
### http://savannah.gnu.org/bugs/?func=detailbug&bug_id=1517&group_id=71

default:   all

MAKEFILE/ :=$(dir $(firstword $(MAKEFILE_LIST)))

ifneq ($(shell ls -d ${MAKEFILE/}/Make.compilers 2> /dev/null), ${MAKEFILE/}/Make.compilers)
  MAKECOMPILERSNOTHING := $(shell ${MAKEFILE/}configure.sh)
endif

include ${MAKEFILE/}Make.compilers


##### non-recursive make magic
# all directories which have been included
//           :=
# current path (empty)
/            :=

# recursive directories to be filtered out of //
# and handled differently
//-RECURSIVE    :=
define MakeRecursive
//-RECURSIVE :=$$/.
endef

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
ifndef WITHOUT
  WITHOUT:=
endif
WITHOUT_:=$(patsubst %,%.,$(strip ${WITHOUT}))

## First Make.include inclusion
# invoke the toplevel include file.
# We use 'Include_File' instead of 'Include' since $/ is empty
$(eval $(call Include_File,$/))

#### Targets which have been declared RECURSIVE are removed
# from // and processed separately.
//            :=$(filter-out ${//-RECURSIVE},${//})

# //-RECURSIVE now holds the paths which are being done legacy style
# and // holds the paths which are going to be part of the system.

### Building subgoals in Make.rules
# subgoals and depends are done with the deferred '=' not the ':='.
# This is because we want dynamic scoping.
__SUBGOALS__=
__DEPGOALS__=

-include ${MAKEFILE/}Make.rules

# now we bring in the depends files as defined by the Includes
# and the patterns in Make.rules
$(eval DEPENDS:=$(foreach x,${//},$(call __DEPGOALS__,$x)))
ifneq ($(strip ${DEPENDS}),)
  # this conditional gets us a way out if things go way wrong
  ifeq ($(filter %-clean,${MAKECMDGOALS}),)
    -include ${DEPENDS} 
  endif
endif

### the standard make targets, applied to all subdirectory targets
# We define the basic form of the all, clean, ... rules on a
# per-path basis (the $/.all, $..clean, ... targets).  This allows
# selective targeting.
## rules for each subtarget
#  current subtargets are :
#        all, ls, clean, real-clean, depends-clean, install
#  with all, clean, real-clean, install being required targets
#  of recursive makes.

clean:         ${//-RECURSIVE:.=.clean}      $/.clean
define .RULE-clean
${1:.=.clean}: $${${1:.=.SUBS}:.=.clean}
	${RM} $${${1:.=.CLEAN}} ${__SUBGOALS__}
ifneq ($(strip ${C_TMP_COMPILE} ${CXX_TMP_COMPILE}),)
	(cd $1 && ${RM} -r ${C_TMP_COMPILE} ${CXX_TMP_COMPILE})
endif

endef
$(eval $(foreach x,${//},$(call .RULE-clean,$x)))


depends-clean:                            $/.depends-clean
${//-RECURSIVE:.=.depends-clean}:
define .RULE-depends-clean
${1:.=.depends-clean}: $${${1:.=.SUBS}:.=.depends-clean}
	${RM} ${1:.=Make.depends} ${__DEPGOALS__}

endef
$(eval $(foreach x,${//},$(call .RULE-depends-clean,$x)))

real-clean:    ${//-RECURSIVE:.=.real-clean} $/.real-clean
define .RULE-real-clean
${1:.=.real-clean}: $${${1:.=.SUBS}:.=.real-clean}
	${RM} $${${1:.=.CLEAN}} ${__SUBGOALS__} ${1:.=Make.depends} ${__DEPGOALS__} $${${1:.=.REAL-CLEAN}}
#ifneq ($(strip ${INSTALL/}),)
#	${RM} -r ${INSTALL/}
#endif
ifneq ($(strip ${C_TMP_COMPILE} ${CXX_TMP_COMPILE}),)
	(cd $1 && ${RM} -r ${C_TMP_COMPILE} ${CXX_TMP_COMPILE})
endif

endef
$(eval $(foreach x,${//},$(call .RULE-real-clean,$x)))

all:           ${//-RECURSIVE:.=.all}        $/.all
define .RULE-all
${1:.=.all}: $${${1:.=.SUBS}:.=.all} ${__SUBGOALS__}

endef
$(eval $(foreach x,${//},$(call .RULE-all,$x)))

## INSTALL:
# Define targets which give a basic recursive traversal to hang
# things from for anything install related that I haven't thought of

# this does the copy part of the install
${//-RECURSIVE:.=.install-copy}:
install-copy:       ${//-RECURSIVE:.=.install-copy}    $/.install-copy
define .RULE-install-copy
${1:.=.install-copy}: $${${1:.=.SUBS}:.=.install-copy}

endef
$(eval $(foreach x,${//},$(call .RULE-install-copy,$x)))

# define the actual install target as being a combination of the
# all target plus the pre-/install-copy targets

install:       ${//-RECURSIVE:.=.install}    $/.install
define .RULE-install
${1:.=.install}: ${1:.=.all} ${1:.=.install-copy}

endef
$(eval $(foreach x,${//},$(call .RULE-install,$x)))



# a sort of debugging thing.  Let's me check on which targets
# actually didn't get made, or got partially made
LSOPTS:=-l
ls:                                       $/.ls
${//-RECURSIVE:.=.ls}:
define .RULE-ls
${1:.=.ls}: $${${1:.=.SUBS}:.=.ls}
	@ files='$$(strip ${__SUBGOALS__})'; \
	  if [ -n "$$$${files}" ] ; then \
	    ls ${LSOPTS} $$$${files} ; \
	  fi ; exit 0;

endef
$(eval $(foreach x,${//},$(call .RULE-ls,$x)))




