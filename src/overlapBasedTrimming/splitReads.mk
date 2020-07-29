TARGET   := splitReads
SOURCES  := splitReads.C \
            splitReads-workUnit.C \
            splitReads-subReads.C \
            splitReads-trimBad.C \
            adjustNormal.C \
            adjustFlipped.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a
