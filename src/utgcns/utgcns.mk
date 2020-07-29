TARGET   := utgcns
SOURCES  := utgcns.C stashContains.C unitigConsensus.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../overlapInCore/libedlib libpbutgcns libboost

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
