# -*- makefile -*-

build:=$(shell while [ `pwd` != "/" -a ! -d build ] ; do \
                 cd .. ; \
               done ; \
               echo `pwd` )

ifeq (${build},/)
  $(error Unable to find the build directory cds/IR/build. Check out this directory.)
endif

build:=${build}/build/

all clean:
	@ echo "(Build directory found at ${build})"
	@ echo "* If you are in a hurry:"
	@ echo "** If using bash then"
	@ echo "**    alias make='/usr/local/ir/bin/gmake -f ${build}Makefile'"
	@ echo "** If using tcsh then"
	@ echo "**    alias make '/usr/local/ir/bin/gmake -f ${build}Makefile'"
	@ echo "** and then do"
	@ echo "**    make"
	@ echo "* If you have some time:"
	@ echo "**    cd ${build}"
	@ echo "**    cvs update"
	@ echo "**    more NOTES"

