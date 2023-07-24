TARGET   := utgcns
SOURCES  := utgcns.C \
            stashContains.C \
            unitigConsensus.C \
            unitigPartition.C \
            utgcns-parameters.C \
            utgcns-processTigs.C \
            libpbutgcns/AlnGraphBoost.C

SRC_INCDIRS := ../utility/src ../stores libpbutgcns libboost

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
