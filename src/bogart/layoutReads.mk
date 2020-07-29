TARGET   := layoutReads
SOURCES  := layoutReads.C \
            AS_BAT_Instrumentation.C \
            AS_BAT_Logging.C \
            AS_BAT_OptimizePositions.C \
            AS_BAT_Outputs.C \
            AS_BAT_OverlapCache.C \
            AS_BAT_PlaceContains.C \
            AS_BAT_PlaceReadUsingOverlaps.C \
            AS_BAT_ReadInfo.C \
            AS_BAT_SetParentAndHang.C \
            AS_BAT_TigVector.C \
            AS_BAT_Unitig.C \
            AS_BAT_Unitig_AddRead.C \
            AS_BAT_Unitig_PlaceReadUsingEdges.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
