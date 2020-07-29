TARGET   := splitHaplotype
SOURCES  := splitHaplotype.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../utgcns

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a
