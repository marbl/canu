TARGET   := tgStoreCompress
SOURCES  := tgStoreCompress.C

SRC_INCDIRS := .. ../utility/src/utility

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a
