TARGET   := errorEstimate
SOURCES  := errorEstimate.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../overlapInCore ../utgcns/libNDalign ../overlapErrorAdjustment

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
