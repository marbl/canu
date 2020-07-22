TARGET   := prefixEditDistance-matchLimitGenerate
SOURCES  := prefixEditDistance-matchLimitGenerate.C

SRC_INCDIRS  := ../.. ../../utility/src/utility ../../stores

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
