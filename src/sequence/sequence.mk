
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)
endif

#            sequence-shiftregister-emit-fast.C \
#            sequence-shiftregister-emit-slow.C \

TARGET   := sequence
SOURCES  := sequence.C \
            sequence-extract.C \
            sequence-generate.C \
            sequence-mutate.C \
            sequence-sample.C \
            sequence-shiftregister.C \
            sequence-shiftregister-emit-fast.C \
            sequence-shiftregister-search-fast.C \
            sequence-shiftregister-search-slow.C \
            sequence-simulate.C \
            sequence-summarize.C

SRC_INCDIRS  := . .. ../utility

#  If we're part of Canu, build with canu support.
#  Otherwise, don't.

ifneq ($(wildcard stores/sqStore.H), )

SRC_CXXFLAGS := -DCANU

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

else

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lmeryl
TGT_PREREQS := libmeryl.a

endif

SUBMAKEFILES :=
