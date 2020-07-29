TARGET   := alignGFA
SOURCES  := alignGFA.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../overlapInCore/libedlib

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
