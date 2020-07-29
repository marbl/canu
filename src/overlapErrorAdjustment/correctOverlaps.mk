TARGET   := correctOverlaps
SOURCES  := correctOverlaps.C \
            correctOverlaps-Correct_Frags.C \
            correctOverlaps-Read_Olaps.C \
            correctOverlaps-Redo_Olaps.C \
            correctOverlaps-Prefix_Edit_Distance.C

SRC_INCDIRS  := .. ../utility/src/utility ../stores ../overlapInCore/liboverlap

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
