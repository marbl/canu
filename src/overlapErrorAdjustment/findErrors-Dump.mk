TARGET   := findErrors-Dump
SOURCES  := findErrors-Dump.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../overlapInCore/liboverlap

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
