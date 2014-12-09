#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../../$(OSTYPE)-$(MACHINETYPE)/bin
endif

#  This will eventually hold kmer's libkmer and libmeryl

#TARGET   := libCA.a

#SOURCES  := 

#SRC_INCDIRS  := .. ../AS_UTL

SUBMAKEFILES := libleaff.mk \
                leaff.mk \
                meryl.mk \
                estimate-mer-threshold.mk

#                mercy.mk \
