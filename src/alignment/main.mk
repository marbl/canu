#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := libCA.a

SOURCES  := alignment-drivers.C \
            analyze-affine-overlap.C \
            brute-force-dp.C \
            dp-align.C \
            find-local-overlap.C \
            find-local-segments.C \
            fix-overlapping-pieces.C \
            local-overlapper.C \
            overlap.C

SRC_INCDIRS  := .. ../AS_UTL ../stores
SUBMAKEFILES := 
