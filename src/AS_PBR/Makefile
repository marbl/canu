#   Copyright (C) 2011, Battelle National Biodefense Institute (BNBI);
#   all rights reserved. Authored by: Sergey Koren
#   
#   This Software was prepared for the Department of Homeland Security
#   (DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
#   part of contract HSHQDC-07-C-00020 to manage and operate the National
#   Biodefense Analysis and Countermeasures Center (NBACC), a Federally
#   Funded Research and Development Center.
#   
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#   
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   
#   * Neither the name of the Battelle National Biodefense Institute nor
#     the names of its contributors may be used to endorse or promote
#     products derived from this software without specific prior written
#     permission.
#   
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
###########################################################################
#

LOCAL_WORK = $(shell cd ../..; pwd)

#  SAMTOOLS should point to a configured, compiled and installed copy of
#  samtools
#
#  We attempt to find it in the same place that our src/ directory is.


#  We'd like to include the global include (c_make.as) to get these
#  definitions, but we cannot without defining rules.

OSTYPE      = $(shell echo `uname`)
MACHINETYPE = $(shell echo `uname -m`)

ifeq ($(MACHINETYPE), x86_64)
  MACHINETYPE = amd64
endif
ifeq ($(MACHINETYPE), Power Macintosh)
  MACHINETYPE = ppc
endif

ifneq "$(origin SAMTOOLS)" "environment"
   ifeq ($(shell ls -d $(LOCAL_WORK)/samtools 2> /dev/null), $(LOCAL_WORK)/samtools)
      SAMTOOLS = $(LOCAL_WORK)/samtools/
   endif
endif

PBR_SOURCES = OutputLayout.cc TrimFastqByQVWindow.cc CorrectPacBio.cc AS_PBR_dumpStore.cc
ifdef SAMTOOLS
   PBR_SOURCES += ConvertSamToCA.cc 
   INC_IMPORT_DIRS += $(SAMTOOLS)
endif

PBR_OBJECTS = $(PBR_SOURCES:.cc=.o)

PBR_LIB_SOURCES = \
				AS_PBR_util.cc \
				AS_PBR_store.cc \
				AS_PBR_correct.cc \
				AS_PBR_filter.cc \
				AS_PBR_mates.cc \
				AS_PBR_output.cc
PBR_LIB_OBJECTS = $(PBR_LIB_SOURCES:.cc=.o)

SOURCES     = $(PBR_SOURCES) $(PBR_LIB_SOURCES)
OBJECTS     = $(PBR_OBJECTS) $(PBR_LIB_OBJECTS)

LIBRARIES     = libAS_PBR.a libCA.a
LIBS          = libCA.a
SAM_LIBS      = $(LIBS) libbam.a

CXX_PROGS = outputLayout trimFastqByQVWindow correctPacBio dumpPBRLayoutStore
ifdef SAMTOOLS
   CXX_PROGS += convertSamToCA
endif

SCRIPTS     = PBcR convertToPBCNS addCNSToStore mhap falcon_sense pbutgcns jellyfish

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

ifeq ($(OSTYPE), Linux)
  LDFLAGS  += -lpthread -lz
  CFLAGS   += -pthread
  CPPFLAGS += -pthread
endif

ifeq ($(OSTYPE), OSF1)
  LDFLAGS  += -lpthread -lz 
  CFLAGS   += -pthread
  CPPFLAGS += -pthread
endif

ifeq ($(OSTYPE), FreeBSD)
   LDFLAGS += -lz
endif

ifeq ($(OSTYPE), Darwin)
   LDFLAGS += -lz
endif


all:    $(OBJECTS) $(LIBRARIES) $(CXX_PROGS) $(SCRIPTS)
	@test -n nop

libAS_PBR.a: $(PBR_LIB_OBJECTS)

libCA.a: $(PBR_LIB_OBJECTS)

libbam.a:
	$(MAKE) -C $(SAMTOOLS)/
	@cp $(SAMTOOLS)/$@ $(LOCAL_LIB)/$@


outputLayout:    OutputLayout.o $(LIBS)

trimFastqByQVWindow:	TrimFastqByQVWindow.o $(LIBS)

convertSamToCA:		ConvertSamToCA.o $(SAM_LIBS)

correctPacBio:	   CorrectPacBio.o $(LIBS)

dumpPBRLayoutStore:	AS_PBR_dumpStore.o $(LIBS)

PBcR: PBcR.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@
	@rm -f $(LOCAL_BIN)/pacBioToCA
	@ln -s $(LOCAL_BIN)/$@ $(LOCAL_BIN)/pacBioToCA

convertToPBCNS: convertToPBCNS.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

addCNSToStore: addCNSToStore.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

# for now we support these as binaries, we should do better and build code

#  This is _ALWAYS_ unpacked, regardless of up-to-dateness.
mhap: mhap.tar.gz
	@tar xvzf $< -C $(LOCAL_LIB)/

falcon_sense.$(LOCAL_NAME).bin:
	@touch falcon_sense.$(LOCAL_NAME).bin

falcon_sense: falcon_sense.$(LOCAL_NAME).bin
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

pbutgcns.$(LOCAL_NAME).bin:
	@touch pbutgcns.$(LOCAL_NAME).bin

pbutgcns: pbutgcns.$(LOCAL_NAME).bin
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

jellyfish.$(LOCAL_NAME).bin:
	@touch jellyfish.$(LOCAL_NAME).bin

ifeq ($(OSTYPE), Darwin)
jellyfish: jellyfish.$(LOCAL_NAME).bin
	@cp $< $(LOCAL_BIN)/$@
	@cp *dylib $(LOCAL_LIB)/
	@chmod 775 $(LOCAL_BIN)/$@
else
jellyfish: jellyfish.$(LOCAL_NAME).bin
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@
endif
