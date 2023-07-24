TARGET   := splitReads
SOURCES  := splitReads.C \
            splitReads-workUnit.C \
            splitReads-subReads.C \
            splitReads-trimBad.C \
            adjustNormal.C \
            adjustFlipped.C

SRC_INCDIRS  := ../utility/src ../stores

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
