/*
Copyright (C) 2011, Battelle National Biodefense Institute (BNBI);
all rights reserved. Authored by: Sergey Koren

This Software was prepared for the Department of Homeland Security
(DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
part of contract HSHQDC-07-C-00020 to manage and operate the National
Biodefense Analysis and Countermeasures Center (NBACC), a Federally
Funded Research and Development Center.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of the Battelle National Biodefense Institute nor
  the names of its contributors may be used to endorse or promote
  products derived from this software without specific prior written
  permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef AS_PBR_STORE_H
#define AS_PBR_STORE_H

#include "AS_PBR_util.hh"
#include "AS_PER_genericStore.h"

#include <map>
#include "boost/dynamic_bitset.hpp"

static const char *rcsid_AS_PBR_STORE_H = "$Id: AS_PBR_store.hh,v 1.1 2012-08-20 13:10:37 skoren Exp $";

/*
 * Store and structures to read/write short-read mapping information
 *
 * This store is responsible for efficiently answering the question, which PacBio reads does a short read map to
 * Using a bit structure, it can also answer the question, what PacBio sequences do these short-reads share
 */
class ShortMapRecord {
public:
	ShortMapRecord(AS_IID offset, AS_IID numRecords, uint16 mapMappings);
	~ShortMapRecord();
	AS_IID readIID;

	bool existsID(AS_IID &id) {
		return mappedBitSet->test(id - bitSetOffset);
	}
	void addMapping(OverlapPos &olap);
	void clearMappings();
	void initMappings(uint32 maxMappings);

	boost::dynamic_bitset<> * getMappedReads() {
		return mappedBitSet;
	}
	SeqInterval *getMapping(uint32 i) {
		if (mappedBitSet->test(i) == false) {
			return NULL;
		}
		uint16 pos = getMappedPosition(i);
		assert(pos >= 0 && pos < mappingSize);
		return &(mappings[pos]);
	}

private:
	ShortMapRecord();
	SeqInterval *mappings;
	//OverlapPos *mappedUnsorted;
	map<AS_IID, uint16> *positionMap;
	uint16 mappingSize;
	bool	isSorted;

	boost::dynamic_bitset<> *mappedBitSet;
	uint32 bitSetBits;
	uint32 bitSetOffset;

	uint64 varStoreOffset;

	inline uint16 getMappedPosition(uint32 i) {
		assert(mappedBitSet != NULL);
		assert(mappings != NULL);
		map<AS_IID, uint16>::const_iterator iter = positionMap->find(i);
		assert(iter != positionMap->end());

		return iter->second;
	}

	friend class ShortMapStore;
};

class ShortMapStore {
public:
	ShortMapStore(const char *path, bool writable, bool creating, bool mem = false);
	~ShortMapStore();

	void appendRecord(ShortMapRecord &mapRecord);
	ShortMapRecord* getRecord(AS_IID index);
	void unlink();

	AS_IID getMappedIID(uint32 i) {
		assert(bitSetOffset != 0);
		return i + bitSetOffset;
	}

	uint32 getStoreIID(AS_IID i) {
		assert(bitSetOffset != 0);
		return i - bitSetOffset;
	}

private:
	ShortMapStore();
	void flush(uint64 neededBytes);

	char                     storePath[FILENAME_MAX];
    StoreStruct             *mapStore;  //  fixed-length fragment info
    FILE					*varFile;   // variable-length info
    bool					deleted;
    HashTable_AS            *IIDtoindex;
    bool					lastWasWrite;
	uint32 					bitSetOffset;

	// buffering data structures
	bool					inMemory;
	uint64					varFileOffset;
	uint64					varBufferLen;
	uint64					varBufferPos;
	uint64					varBufferMax;
	char 					*varBuffer;
};

/*
 * Store and structures to read/write intermediate layout format
 *
 * This store is responsible for storing information on a tiling for each PacBio sequence
 */
struct LayRecordStore {
  int           bufferLen;    //  length of valid data in the buffer
  int           bufferPos;    //  position the read is at in the buffer
  int           bufferMax;    //  allocated size of the buffer
  bool			isOutput;		// whether we are reading or writing
  bool			isPopened;		// whether we are reading or writing from a pipe and not a file
  int32       *buffer;
  FILE         *file;
};

static const uint32 LAY_RECORD_SIZE = 1 + 4 + AS_OVS_NWORDS;
static const uint32 LAY_HEADER_SIZE = 2;

struct LayRecord {
	AS_IID iid;							// the iid of the PacBio sequence
	vector<OverlapPos> mp;				// the tiling information for each short-read
    map<AS_IID, SeqInterval> bClrs;		// the clear-ranges for each short-read
    map<AS_IID, OVSoverlap> bOvls;		// the overlaps that generated this tiling
    map<AS_IID, int32> readStarts;		// the starting position for each sequence in the tiling (optionally populated)
    map<AS_IID, bool> readOri;			// the orientations of each sequence in the tiling (optionally populated)
    uint32 length;						// the length of the PacBio sequence (computed based on the tiling)
};

static void clearLayout(LayRecord &r) {
	r.mp.clear();
	r.bClrs.clear();
	r.bOvls.clear();
	r.readStarts.clear();
	r.readOri.clear();
	r.length = 0;
	r.iid = 0;
}

extern LayRecordStore *createLayFile(const char *name);
extern LayRecordStore *openLayFile(const char *name);
extern void	closeLayFile(LayRecordStore *lof);
extern uint32 writeLayRecord(LayRecordStore *out, LayRecord &r, boost::dynamic_bitset<> *bits, double percentageToStoreInBits, map<AS_IID, bool> *subset = NULL, BinaryOverlapFile *bof = NULL);
extern bool readLayRecord(LayRecordStore *in, LayRecord &r);
extern bool readLayRecord(LayRecordStore *in, LayRecord &r, PBRThreadGlobals *waGlobal);
extern bool readLayRecord(LayRecordStore *in, LayRecord &r, MultiAlignT *ma);
#endif
