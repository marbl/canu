TARGET   := sqStoreDumpFASTQ
SOURCES  := sqStoreDumpFASTQ.C

SRC_INCDIRS := .. ../stores ../utility/src/utility ../overlapBasedTrimming

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
