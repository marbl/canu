# -*- makefile -*-

### NOTE: make-3.80 has a bug in eval which causes it not to return
###   an empty target, but instead it sometimes returns very weird
###   stuff that makes make segv, I think.
### You need to get it patched at:
### http://savannah.gnu.org/bugs/?func=detailbug&bug_id=1517&group_id=71

default:   all

# kind of useful
debug:
	@echo DIRS ${//}
	@echo RECS ${RECURSIVE}
	@echo C_EXES ${C_EXES}
	@echo C_LIBS ${C_LIBS}
	@echo CXX_EXES ${CXX_EXES}
	@echo CXX_LIBS ${CXX_LIBS}
	@echo CXX_SHLIBS ${CXX_SHLIBS}
	@echo INSTALL ${INSTALL/}

LSOPTS:=-l
ls:
	@ls ${LSOPTS} ${TARGETS}

# Compiler information and options.  Always get it from the same directory
# as this makefile.
MAKEFILE/ :=$(dir $(firstword $(MAKEFILE_LIST)))
-include ${MAKEFILE/}Make.compilers

##### non-recursive make magic
//           :=
RECURSIVE    :=
/            :=

define Include
  ifeq ($(wildcard ${1}Make.include),${1}Make.include)
    $/.SUBS +=${1}.
    // +=${1}.
    ${1}.SUBS :=
    /  :=${1}
    include ${1}Make.include
    /  :=$/
  endif
endef

##### System specific includes

# invoke the toplevel include file
$(eval $(call Include,$/))

#### Pulling things up the make tree ####

//            :=$(filter-out ${RECURSIVE},${//})

C_LIBS       :=${foreach d,${//},${${d:.=.C_LIBS}}}
C_EXES       :=${foreach d,${//},${${d:.=.C_EXES}}}
CXX_LIBS     :=${foreach d,${//},${${d:.=.CXX_LIBS}}}
CXX_EXES     :=${foreach d,${//},${${d:.=.CXX_EXES}}}
CXX_SHLIBS   :=${foreach d,${//},${${d:.=.CXX_SHLIBS}}}
TARGETS      :=${C_LIBS} ${C_EXES} ${CXX_LIBS} ${CXX_EXES} ${CXX_SHLIBS}

### this gives us a way out if something goes horribly wrong...
ifeq ($(filter %real-clean,${MAKECMDGOALS}),)
  include     ${//:.=Make.depends}
endif
## rules for each subtarget
#  current subtargets are : all, clean, real-clean, depends-clean

define .RULE-clean
${1:.=.clean}: $${${1:.=.SUBS}:.=.clean}
	rm -f \
	  $${${1:.=.CLEAN}} \
	  $${${1:.=.C_LIBS}} $${${1:.=.C_EXES}} \
	  $${${1:.=.CXX_LIBS}} $${${1:.=.CXX_EXES}} $${${1:.=.CXX_SHLIBS}}
	(cd $1 && \
         rm -rf ${C_TMP} ${C_TMP_ARCH} ${CXX_TMP} ${CXX_TMP_ARCH})

endef
define .RULE-depends-clean
${1:.=.depends-clean}: $${${1:.=.SUBS}:.=.depends-clean}
	rm -f ${1:.=Make.depends}

endef
define .RULE-real-clean
${1:.=.real-clean}: ${1:.=.clean} ${1:.=.depends-clean}
	rm -f $${${1:.=.REAL-CLEAN}}

endef
define .RULE-all
${1:.=.all}: $${${1:.=.SUBS}:.=.all} \
   $${${1:.=.C_EXES}} $${${1:.=.CXX_EXES}} \
   $${${1:.=.C_LIBS}} $${${1:.=.CXX_LIBS}} \
   $${${1:.=.CXX_SHLIBS}}

endef
define .RULE-install
${1:.=.install}: $${${1:.=.SUBS}:.=.install} \
   $${${1:.=.C_EXES}} $${${1:.=.CXX_EXES}} \
   $${${1:.=.C_LIBS}} $${${1:.=.CXX_LIBS}} \
   $${${1:.=.CXX_SHLIBS}} \
   $${${1:.=.C_INCS}} $${${1:.=.CXX_INCS}}
	@echo "***for $$@"
	List='$$(strip $${${1:.=.C_EXES}} $${${1:.=.CXX_EXES}})'; \
	if [ -n "$$$${List}" ] ; then \
	  mkdir -p ${INSTALL/}bin ; \
	  cp -p $$$${List} ${INSTALL/}bin ; \
        fi
	List='$$(strip $${${1:.=.C_LIBS}} $${${1:.=.CXX_LIBS}} $${${1:.=.CXX_SHLIBS}})'; \
	if [ -n "$$$${List}" ] ; then \
	  mkdir -p ${INSTALL/}lib/$${${1:.=.LIB_SUB/}} ; \
	  cp -p $$$${List} ${INSTALL/}lib/$${${1:.=.LIB_SUB/}} ; \
        fi
	List='$$(strip $${${1:.=.C_INCS}} $${${1:.=.CXX_INCS}})'; \
	if [ -n "$$$${List}" ] ; then \
	  mkdir -p ${INSTALL/}include/$${${1:.=.INC_SUB/}} ; \
	  cp -p $$$${List} ${INSTALL/}include/$${${1:.=.INC_SUB/}} ; \
        fi

endef

# each rule applies within its local subdirectory
$(eval $(foreach x,${//},$(call .RULE-all,$x)))
$(eval $(foreach x,${//},$(call .RULE-clean,$x)))
$(eval $(foreach x,${//},$(call .RULE-depends-clean,$x)))
$(eval $(foreach x,${//},$(call .RULE-real-clean,$x)))
$(eval $(foreach x,${//},$(call .RULE-install,$x)))

## some general targets must be supported for legacy recursive-make
#  these targets will come before the non-recursive targets.
#  Maybe it would be better to have before and after recmakedirs.
recursive-all:           ${RECURSIVE:.=.all}
recursive-clean:         ${RECURSIVE:.=.clean}
recursive-real-clean:    ${RECURSIVE:.=.real-clean}
recursive-install:       ${RECURSIVE:.=.install}
# depends clean is not a recursive target
${RECURSIVE:.=.depends-clean}:

###
# finally, the global make targets
all:           recursive-all           $/.all
clean:         recursive-clean         $/.clean
real-clean:    recursive-real-clean    $/.real-clean
depends-clean:                         $/.depends-clean
install:       recursive-install       $/.install

##  generic pattern rules
%.o: %.c
	${CC} ${CFLAGS} ${CFLAGS_ARCH} -o $@ -c $<

%.o: %.cc
	${CXX} ${CXXFLAGS} ${CXXFLAGS_ARCH} -o $@ -c $<

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${CXXFLAGS_ARCH} -o $@ -c $<

%.o: %.C
	${CXX} ${CXXFLAGS} ${CXXFLAGS_ARCH} -o $@ -c $<


## EXE targets
${C_EXES}:
	${CC} ${CLDFLAGS} ${CLDFLAGS_ARCH} -o $@ $^ ${CLIBS} ${CLIBS_ARCH}

${CXX_EXES}:
	${CXX} ${CXXLDFLAGS} ${CXXLDFLAGS_ARCH} -o $@ $^ ${CXXLIBS} ${CXXLIBS_ARCH}

## LIB targets
${C_LIBS}:
	rm -f $@
	${AR} ${ARFLAGS} $@ $^

${CXX_LIBS}:
	rm -f $@
	${AR} ${ARFLAGS} $@ $^

## Shared targets
${CXX_SHLIBS}:
	${CXX} ${CXXFLAGS} ${CXXFLAGS_ARCH} ${SHLIB_FLAGS} -o $@ $^ ${CXXLIBS} ${CXXLIBS_ARCH}


## Make.depends dependency calculation
${//:.=Make.depends}:
	@ echo "(creating $@)"
	@ echo "# -*- makefile -*-" > $@
	@ echo "# Make.depends file: $@" >> $@
	\
	List='$(strip ${${@:Make.depends=.C_SRCS}})'; \
	if [ -n "$${List}" ] ; then \
	  for f in $${List} ; do \
	    dir=`dirname $${f}`/ ; \
	    dir=`echo $${dir} | sed -e's!^\./!!'`; \
	    ${CCDEP} ${CFLAGS} $${f} \
	    | sed -e"/:/s!^!$${dir}!" \
	    >> $@ ; \
	  done \
	fi
	\
	List='$(strip ${${@:Make.depends=.CXX_SRCS}})'; \
	if [ -n "$${List}" ] ; then \
	  for f in $${List} ; do \
	    dir=`dirname $${f}`/ ; \
	    dir=`echo $${dir} | sed -e's!^\./!!'`; \
	    ${CXXDEP} ${CXXFLAGS} $${f} \
	    | sed -e"/:/s!^!$${dir}!" \
	    >> $@ ; \
	  done \
	fi
