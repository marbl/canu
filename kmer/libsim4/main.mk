
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

#CXX = g++46
#CC  = gcc46
#CXXFLAGS := -fopenmp -D_GLIBCXX_PARALLEL -O3 -fPIC -m64 -pipe -Wno-write-strings
#LDFLAGS  := -fopenmp -lm

TARGET   := libsim4.a

SOURCES  := sim4core/sim4command.C \
            sim4core/sim4parameters.C \
            sim4core/sim4string.C \
            sim4core/Xtend1.C \
            sim4core/align.C \
            sim4core/exon_cores.C \
            sim4core/extend.C \
            sim4core/glimmerSplice.C \
            sim4core/greedy.C \
            sim4core/mspManager.C \
            sim4core/pluri_align.C \
            sim4core/poly.C \
            sim4core/sim4b1.C \
            sim4core/sim4b1a.C \
            sim4core/sim4b1-1.C \
            sim4core/sim4b1-2.C \
            sim4core/sim4b1-3.C \
            sim4core/sim4b1-4.C \
            sim4core/sim4b1_s.C \
            sim4core/sites.C \
            sim4core/sites_donor.C \
            sim4core/sites_acceptor.C \
            sim4core/sites_score.C \
            sim4core/splice.C \
            sim4core/table.C \
            sim4core/util.C \
            sim4polish/sim4polish-compare.C \
            sim4polish/sim4polish-copy.C \
            sim4polish/sim4polish-deleteexon.C \
            sim4polish/sim4polish-exons.C \
            sim4polish/sim4polish-polishtostring.C \
            sim4polish/sim4polish-read.C \
            sim4polish/sim4polish-stringtopolish.C \
            sim4polish/sim4polish-updatescores.C \
            sim4polish/sim4polish.C \
            sim4polish/sim4polishList.C \
            sim4polish/sim4polishBuilder.C \
            sim4polish/sim4polishFile.C \
            sim4polish/sim4polishReader.C \
            sim4polish/sim4polishWriter.C

SRC_INCDIRS  := ../libutil ../libbio ../libseq sim4core sim4polish

SUBMAKEFILES := 
