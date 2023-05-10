TARGET   := loadTrimmedReads
SOURCES  := loadTrimmedReads.C

SRC_INCDIRS := .. ../stores ../utility/src/utility ../overlapBasedTrimming

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
