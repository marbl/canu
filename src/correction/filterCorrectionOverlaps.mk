TARGET   := filterCorrectionOverlaps
SOURCES  := filterCorrectionOverlaps.C

SRC_INCDIRS  := .. ../utility/src ../stores

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
