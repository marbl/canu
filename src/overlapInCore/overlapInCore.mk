TARGET   := overlapInCore
SOURCES  := overlapInCore.C \
            overlapInCore-Build_Hash_Index.C \
            overlapInCore-Find_Overlaps.C \
            overlapInCore-Output.C \
            overlapInCore-Process_Overlaps.C \
            overlapInCore-Process_String_Overlaps.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores liboverlap

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
