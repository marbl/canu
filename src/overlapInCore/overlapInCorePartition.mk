TARGET   := overlapInCorePartition
SOURCES  := overlapInCorePartition.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores liboverlap

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
