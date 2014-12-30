#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := unitigger
SOURCES  := AS_CGB_Bubble.C \
            AS_CGB_Bubble_Graph.C \
            AS_CGB_Bubble_Popper.C \
            AS_CGB_Bubble_PopperMethods.C \
            AS_CGB_Bubble_VertexSet.C \
            AS_CGB_Bubble_dfs.C \
            AS_CGB_cga.C \
            AS_CGB_cgb.C \
            AS_CGB_chimeras.C \
            AS_CGB_classify.C \
            AS_CGB_edgemate.C \
            AS_CGB_fga.C \
            AS_CGB_fgb.C \
            AS_CGB_histo.C \
            AS_CGB_main.C \
            AS_CGB_traversal.C \
            AS_CGB_unitigger.C \
            AS_CGB_walk.C \
            AS_FGB_contained.C \
            AS_FGB_hanging_fragment.C \
            AS_FGB_io.C \
            AS_FGB_main.C

SRC_INCDIRS  := .. ../AS_UTL ../stores ../alignment

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lCA
TGT_PREREQS := libCA.a

SUBMAKEFILES :=
