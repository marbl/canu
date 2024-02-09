TARGET   := utgcns
SOURCES  := utgcns.C \
            stashContains.C \
            unitigConsensus.C \
            unitigPartition.C \
            utgcns-parameters.C \
            utgcns-processTigs.C \
            libpbutgcns/AlnGraphBoost.C


ifeq (${OSTYPE}, Linux)
  BOOST ?= libboost
endif

ifeq (${OSTYPE}, FreeBSD)
  ifeq ($(wildcard /usr/local/include/boost), /usr/local/include/boost)
    BOOST ?= /usr/local/include
  else
    BOOST ?= libboost
  endif
endif

ifeq (${OSTYPE}, Darwin)
  ifeq      ($(wildcard /opt/homebrew/include/boost), /opt/homebrew/include/boost)
    BOOST ?= /opt/homebrew/include
  else ifeq ($(wildcard /usr/local/include/boost),    /usr/local/include/boost)
    BOOST ?= /usr/local/include
  else ifeq ($(wildcard ${BREW}/opt/boost/include),   ${BREW}/opt/boost/include)
    BOOST ?= ${BREW}/opt/boost/include
  else ifeq ($(wildcard ${PORT}/include/boost),       ${PORT}/include/boost)
    #  UNTESTED!
    BOOST ?= ${PORT}/include/boost
  else
    BOOST ?= libboost
  endif
endif

ifeq (${BOOST}, libboost)
  SRC_CXXFLAGS += -DBOOST_NO_AUTO_PTR
endif

SRC_INCDIRS := ../utility/src ../stores libpbutgcns

EXT_INCDIRS := ${BOOST}

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
