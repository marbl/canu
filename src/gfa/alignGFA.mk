TARGET   := alignGFA
SOURCES  := alignGFA.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../overlapInCore/libedlib

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
