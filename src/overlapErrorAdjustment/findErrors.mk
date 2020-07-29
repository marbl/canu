TARGET   := findErrors
SOURCES  := findErrors.C \
            findErrors-Analyze_Alignment.C \
            findErrors-Output.C \
            findErrors-Prefix_Edit_Distance.C \
            findErrors-Process_Olap.C \
            findErrors-Read_Frags.C \
            findErrors-Read_Olaps.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../overlapInCore/liboverlap

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
