TARGET   := sqStoreDumpMetaData
SOURCES  := sqStoreDumpMetaData.C ../seqrequester/src/seqrequester/summarize.C

SRC_INCDIRS := ../utility/src ../seqrequester/src/seqrequester

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
