TARGET   := errorEstimate
SOURCES  := errorEstimate.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../overlapInCore ../utgcns/libNDalign ../overlapErrorAdjustment

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
