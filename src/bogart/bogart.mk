#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := bogart
SOURCES  := bogart.C \
            AS_BAT_AssemblyGraph.C \
            AS_BAT_BestOverlapGraph.C \
            AS_BAT_ChunkGraph.C \
            AS_BAT_CreateUnitigs.C \
            AS_BAT_DropDeadEnds.C \
            AS_BAT_Instrumentation.C \
            AS_BAT_Logging.C \
            AS_BAT_MarkRepeatReads.C \
            AS_BAT_MergeOrphans.C \
            AS_BAT_Outputs.C \
            AS_BAT_OverlapCache.C \
            AS_BAT_PlaceContains.C \
            AS_BAT_PlaceReadUsingOverlaps.C \
            AS_BAT_PopulateUnitig.C \
            AS_BAT_PromoteToSingleton.C \
            AS_BAT_ReadInfo.C \
            AS_BAT_SetParentAndHang.C \
            AS_BAT_SplitDiscontinuous.C \
            AS_BAT_TigGraph.C \
            AS_BAT_TigVector.C \
            AS_BAT_Unitig.C \
            AS_BAT_Unitig_AddRead.C \
            AS_BAT_Unitig_PlaceReadUsingEdges.C

SRC_INCDIRS  := .. ../AS_UTL ../stores

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
