TARGET   := bogart
SOURCES  := bogart.C \
            AS_BAT_AssemblyGraph.C \
            AS_BAT_BestOverlapGraph.C \
            AS_BAT_ChunkGraph.C \
            AS_BAT_CreateUnitigs.C \
            AS_BAT_DetectSpurs.C \
            AS_BAT_DropDeadEnds.C \
            AS_BAT_Instrumentation.C \
            AS_BAT_Logging.C \
            AS_BAT_MarkRepeatReads.C \
            AS_BAT_MergeOrphans.C \
            AS_BAT_OptimizePositions.C \
            AS_BAT_Outputs.C \
            AS_BAT_OverlapCache.C \
            AS_BAT_PlaceContains.C \
            AS_BAT_PlaceReadUsingOverlaps.C \
            AS_BAT_PopulateUnitig.C \
            AS_BAT_PromoteToSingleton.C \
            AS_BAT_ReadInfo.C \
            AS_BAT_SetParentAndHang.C \
            AS_BAT_SplitDiscontinuous.C \
            AS_BAT_SplitTig_ByRegions.C \
            AS_BAT_SplitTig_AtReadEnds.C \
            AS_BAT_TigGraph.C \
            AS_BAT_TigVector.C \
            AS_BAT_Unitig.C \
            AS_BAT_Unitig_AddRead.C \
            AS_BAT_Unitig_PlaceReadUsingEdges.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
