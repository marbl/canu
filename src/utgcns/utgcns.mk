TARGET   := utgcns
SOURCES  := utgcns.C stashContains.C unitigPartition.C unitigConsensus.C

SRC_INCDIRS  := .. \
                ../utility/src \
                ../stores \
                libpbutgcns \
                libboost \
                htslib \
                htslib/lib \
                htslib/htscodecs/htscodecs \
                ../merlin \
                ../merlin/build/lib

TGT_LDFLAGS := -L${TARGET_DIR}/lib -L/data/solarsj/tools/verkko/src/canu/src/utgcns/htslib/lib # -Lhtslib/lib # 
TGT_LDLIBS  := -l${MODULE} -lhts
TGT_PREREQS := lib${MODULE}.a
