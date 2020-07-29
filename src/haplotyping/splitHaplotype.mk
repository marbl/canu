TARGET   := splitHaplotype
SOURCES  := splitHaplotype.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../utgcns

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
