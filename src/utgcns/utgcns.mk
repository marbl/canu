TARGET   := utgcns
SOURCES  := utgcns.C stashContains.C unitigConsensus.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../overlapInCore/libedlib libpbutgcns libboost

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
