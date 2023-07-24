TARGET   := trimReads
SOURCES  := trimReads.C \
            trimReads-bestEdge.C \
            trimReads-largestCovered.C

SRC_INCDIRS  := ../utility/src ../stores

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
