TARGET   := tgStoreCompress
SOURCES  := tgStoreCompress.C

SRC_INCDIRS := .. ../utility/src/utility

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
