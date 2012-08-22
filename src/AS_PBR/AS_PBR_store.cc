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

static const char *rcsid_AS_PBR_STORE_C = "$Id: AS_PBR_store.cc,v 1.2 2012-08-22 14:41:00 skoren Exp $";

using namespace std;

#include "AS_PBR_store.hh"
#include "MultiAlignMatePairAnalysis.H"
#include "MultiAlign.h"
#include "AS_UTL_alloc.h"

#include <vector>

static int LOF_BufferMax = 32 * 1024 * 1024;
static int VERBOSE = 0;

ShortMapRecord::ShortMapRecord() {
	// empty record, should only be used for reading from stream
	varStoreOffset = bitSetOffset = 0;
	mappedBitSet = NULL;
	mappings = NULL;
	//mappedUnsorted = NULL;
	isSorted = false;
}

ShortMapRecord::ShortMapRecord(AS_IID offset, AS_IID maxSequenceID, uint16 maxMappings) {
	// empty
	mappings = new SeqInterval[maxMappings + 1];
	//mappedUnsorted = new OverlapPos[maxMappings];
	memset(mappings, 0, sizeof(SeqInterval) * maxMappings);
	//memset(mappedUnsorted, 0, sizeof(OverlapPos) * maxMappings);

	positionMap = new map<AS_IID, uint16>();

	mappingSize = 0;
	varStoreOffset = 0;
	bitSetOffset = offset;

	mappedBitSet = new boost::dynamic_bitset<>(maxSequenceID - bitSetOffset + 1);
	isSorted = false;
}

void
ShortMapRecord::clearMappings() {
	memset(mappings, 0, sizeof(SeqInterval) * mappingSize);
	//memset(mappedUnsorted, 0, sizeof(OverlapPos) * mappingSize);
	mappingSize = 0;
	mappedBitSet->reset();
	positionMap->clear();
}

void
ShortMapRecord::addMapping(OverlapPos &olap) {
	mappedBitSet->set(olap.ident - bitSetOffset, true);
	uint32 pos = mappingSize;
	positionMap->insert(pair<AS_IID, uint16>(olap.ident-bitSetOffset, mappingSize)); //getMappedPosition(olap.ident - bitSetOffset);

	if (!(mappings[pos].bgn == 0 && mappings[pos].end == 0)) {
		// need to insert, move to make room
		memmove(mappings + pos + 1, mappings + pos, sizeof(SeqInterval) * (mappingSize - pos));
	}
	isSorted = true;
	mappings[pos] = olap.position;
	mappingSize++;

	assert(positionMap->size() == mappingSize);
}

ShortMapRecord::~ShortMapRecord() {
	delete[] mappings;
	delete mappedBitSet;
	delete positionMap;
}

ShortMapStore::ShortMapStore() {
	// empty store
}

ShortMapStore::ShortMapStore(const char *path, bool writable, bool creating, bool mem) {
	char  name[FILENAME_MAX];
	char   mode[4];

	bitSetOffset = 0;
	inMemory = mem;

	memset(storePath, 0, FILENAME_MAX);
	strncpy(storePath, path, strlen(path));
	if (writable) {
		strcpy(mode, "r+");
	} else {
		strcpy(mode, "r");
	}

	if (creating) {
		sprintf(name,"%s.shortmap", storePath);
		mapStore = createIndexStore(name, "map", sizeof(ShortMapRecord), 1);
		sprintf(name, "%s.shortmap.var", storePath);
		varFile = fopen(name, "w+");
		sprintf(name, "%s.shortmap.index", storePath);
		IIDtoindex = CreateScalarHashTable_AS();
		SaveHashTable_AS(name, IIDtoindex);
	} else {
		sprintf(name,"%s.shortmap", storePath);
		mapStore   = openStore(name, mode);
		sprintf(name, "%s.shortmap.var", storePath);
		varFile = fopen(name, mode);
		sprintf(name, "%s.shortmap.index", storePath);
		IIDtoindex = LoadUIDtoIIDHashTable_AS(name);
	}
	mapStore = convertStoreToMemoryStore(mapStore);

	// create buffers to avoid writing to disk excessively
	if (LOF_BufferMax % sizeof(SeqInterval) % sizeof(uint16)) {
		fprintf(stderr, "The seq interval must be an even multiple of int32 or we cannot fit it into our store\n");
		assert(0);
	}
	if (LOF_BufferMax % sizeof(SeqInterval) % sizeof(int32)) {
		fprintf(stderr, "The seq interval must be an even multiple of int32 or we cannot fit it into our store\n");
		assert(0);
	}
	if (LOF_BufferMax % sizeof(boost::dynamic_bitset<>::block_type) != 0) {
			fprintf(stderr, "The boost bitset block size must be an even multiple of int32 or we cannot fit it into our store\n");
			assert(0);
	}
	varFileOffset = 0;
	varBufferLen  = 0;
	varBufferPos  = 0;
	varBufferMax  = (LOF_BufferMax);
	varBuffer     = (char*) safe_malloc(varBufferMax);

	deleted = false;
}

void
ShortMapStore::flush(uint64 neededBytes) {
	if (neededBytes == 0 || neededBytes >= varBufferMax - varBufferLen) {
		if (neededBytes != 0 && inMemory) {
			int64 desiredSize = varBufferLen + neededBytes;
			varBufferMax = desiredSize + (LOF_BufferMax);
			varBuffer  = (char *)safe_realloc(varBuffer, varBufferMax);
		} else {
			if (varBufferLen > 0) {
				AS_UTL_safeWrite(varFile, varBuffer, "AS_PBR_flushShortMapStore", sizeof(char), varBufferLen);
				varBufferLen = 0;
			}
			varFileOffset = AS_UTL_ftell(varFile);
		}
	}
}

ShortMapStore::~ShortMapStore() {
	if (!deleted) {
		flush(0);
		closeStore(mapStore);
		fclose(varFile);
		safe_free(varBuffer);

		char  name[FILENAME_MAX];
		sprintf(name, "%s.shortmap.index", storePath);
		SaveHashTable_AS(name, IIDtoindex);
	}
	DeleteHashTable_AS(IIDtoindex);
}

void
ShortMapStore::unlink() {
	char  name[FILENAME_MAX];
	deleted = true;

	delete[] varBuffer;
	closeStore(mapStore);
	fclose(varFile);

	sprintf(name,"%s.shortmap", storePath);
	AS_UTL_unlink(name);
	sprintf(name,"%s.shortmap.var", storePath);
	AS_UTL_unlink(name);
	sprintf(name, "%s.shortmap.index", storePath);
	AS_UTL_unlink(name);
}

void ShortMapStore::appendRecord(ShortMapRecord &mapRecord) {
	// seek to the end if necessary
	if (lastWasWrite == 0) {
		AS_UTL_fseek(varFile, 0, SEEK_END);
	}

	// record information on this record
	mapRecord.varStoreOffset = varFileOffset + varBufferLen;
	mapRecord.bitSetBits = mapRecord.mappedBitSet->size();

	// save the table
	flush((sizeof(SeqInterval) * mapRecord.mappingSize));
	memcpy(varBuffer + varBufferLen, mapRecord.mappings, sizeof(SeqInterval) * mapRecord.mappingSize);
	varBufferLen += sizeof(SeqInterval) * mapRecord.mappingSize;

	// copy the positional map
	assert(mapRecord.positionMap->size() == mapRecord.mappingSize);
	flush((sizeof(AS_IID) + sizeof(uint16)) * mapRecord.mappingSize);
	for (map<AS_IID, uint16>::const_iterator iter = mapRecord.positionMap->begin(); iter != mapRecord.positionMap->end(); iter++) {
		memcpy(varBuffer + varBufferLen, &iter->first, sizeof(iter->first));
		varBufferLen += sizeof(iter->first);
		memcpy(varBuffer + varBufferLen, &iter->second, sizeof(iter->second));
		varBufferLen += sizeof(iter->second);
	}

	appendIndexStore(mapStore, &mapRecord);
	InsertInHashTable_AS(IIDtoindex, mapRecord.readIID, 0, getLastElemStore(mapStore), 0);

	// make sure we correctly stored the number of bits/blocks in the bitset
	assert(ceil((double)mapRecord.bitSetBits / (sizeof(boost::dynamic_bitset<>::block_type) * NUM_BITS)) == mapRecord.mappedBitSet->num_blocks());

	if (VERBOSE) {
	    fprintf (stderr, "The record for %d contains %d maps\n", mapRecord.readIID, mapRecord.mappingSize);
	    for (map<AS_IID, uint16>::const_iterator iter = mapRecord.positionMap->begin(); iter != mapRecord.positionMap->end(); iter++) {
	        fprintf(stderr, "Mapping to %d is from position %d to %d\n", (iter->first + mapRecord.bitSetOffset), mapRecord.mappings[iter->second].bgn, mapRecord.mappings[iter->second].end);
	    }
	}
	lastWasWrite = 1;
}

ShortMapRecord* ShortMapStore::getRecord(AS_IID index) {
	if (ExistsInHashTable_AS(IIDtoindex, index, 0) == 0) {
		return NULL;
	}

	// load the main record
	ShortMapRecord *result = new ShortMapRecord();
	getIndexStore(mapStore, (int32)LookupValueInHashTable_AS(IIDtoindex, index, 0), result);

	// flush if we need to
	if (lastWasWrite) {
		fflush(varFile);
	}

	if (bitSetOffset == 0) {
		bitSetOffset = result->bitSetOffset;
	} else {
		assert(bitSetOffset == result->bitSetOffset);
	}
	// now load variable-length structures
	if (inMemory) {
		if (varBufferLen == 0) {
			AS_UTL_fseek(varFile, 0, SEEK_END);
			uint64 fileSize = AS_UTL_ftell(varFile);
			safe_free(varBuffer);
			varBuffer = (char *) safe_malloc(fileSize);

			varBufferLen = fileSize;
			AS_UTL_fseek(varFile, 0, SEEK_SET);
			AS_UTL_safeRead(varFile, varBuffer, "ShortMapRecord::getRecord", sizeof(int32), varBufferLen);
		}
		uint64 offset = (result->varStoreOffset);
		assert(offset < varBufferLen);
		result->mappings = new SeqInterval[result->mappingSize];
		memcpy(result->mappings, varBuffer + offset, sizeof(SeqInterval) * result->mappingSize);
		offset += sizeof(SeqInterval) * result->mappingSize;

		assert(offset < varBufferLen);
		result->mappedBitSet = new boost::dynamic_bitset<>(result->bitSetBits);

		AS_IID ident;
		uint16 pos;
		result->positionMap = new map<AS_IID, uint16>();
		for (int i = 0; i < result->mappingSize; i++) {
			memcpy(&ident, varBuffer + offset, sizeof(AS_IID));
			offset += sizeof(AS_IID);
			memcpy(&pos, varBuffer + offset, sizeof(uint16));
			offset += sizeof(uint16);
			result->positionMap->insert(pair<AS_IID, uint16>(ident, pos));
			result->mappedBitSet->set(ident, true);
		}
	} else {
		AS_UTL_fseek(varFile, (off_t)result->varStoreOffset, SEEK_SET);
		result->mappings = new SeqInterval[result->mappingSize];
		AS_UTL_safeRead(varFile, result->mappings, "ShortMapStore::getRecord", sizeof(SeqInterval), result->mappingSize);

        // get the map
		result->mappedBitSet = new boost::dynamic_bitset<>(result->bitSetBits);
		result->positionMap = new map<AS_IID, uint16>();
		AS_IID ident;
		uint16 pos;
		for (int i = 0; i < result->mappingSize; i++) {
			AS_UTL_safeRead(varFile, &ident, "ShortMapStore::getRecord", sizeof(AS_IID), 1);
			AS_UTL_safeRead(varFile, &pos, "ShortMapStore::getRecord", sizeof(uint16), 1);

			result->positionMap->insert(pair<AS_IID, uint16>(ident, pos));
            result->mappedBitSet->set(ident, true);
		}

	}
	lastWasWrite = 0;
	return result;
}

static void initializeLOF(LayRecordStore *lof) {
	assert(lof != NULL);

	lof->bufferLen  = 0;
	lof->bufferPos  = 0;
	lof->bufferMax  = ((LOF_BufferMax / (LAY_RECORD_SIZE * sizeof(int32))) * LAY_RECORD_SIZE);
	lof->buffer     = new int32[lof->bufferMax + 1];
	lof->file       = NULL;
	lof->isPopened 	= false;
	lof->isOutput	= false;

	if (lof->bufferMax < (LAY_HEADER_SIZE + LAY_RECORD_SIZE)) {
		fprintf(stderr, "Buffer is too small, can not hold a header with one record at size %d, minimum is %d\n", lof->bufferMax, (LAY_HEADER_SIZE+LAY_RECORD_SIZE));
		assert(0);
	}
}

static void
flushLOF(LayRecordStore *lof) {
	if ((lof->isOutput) && (lof->bufferLen > 0)) {
		AS_UTL_safeWrite(lof->file, lof->buffer, "AS_PBR_flushLOF", sizeof(int32), lof->bufferLen);
		lof->bufferLen = 0;
	}
}

LayRecordStore *
createLayFile(const char *name) {
	char     cmd[1024 + FILENAME_MAX];
	memset(cmd, 0, 1024 + FILENAME_MAX);

	LayRecordStore   *lof = new LayRecordStore;
	initializeLOF(lof);
	lof->isOutput = true;

	if ((name == NULL) || (strcmp(name, "-") == 0))
		name = NULL;

	errno = 0;
	if (name == NULL) {
		lof->file = stdout;
	/*
	 * disable zip file support because it does not seem to work on all filesystem with threads
	} else if (strcasecmp(name+strlen(name)-3, ".gz") == 0) {
		sprintf(cmd, "gzip -1c > %s", name);
		lof->file = popen(cmd, "w");
		lof->isPopened = TRUE;
	} else if (strcasecmp(name+strlen(name)-4, ".bz2") == 0) {
		sprintf(cmd, "bzip2 -1c > %s", name);
		lof->file = popen(cmd, "w");
		lof->isPopened = TRUE;
	*/} else {
		lof->file = fopen(name, "w");
	}
	if (errno) {
		fprintf(stderr, "createLayFile()-- Failed to open '%s' for writing: %s\n",
				name, strerror(errno));
		exit(1);
	}

	return(lof);
}

LayRecordStore *
openLayFile(const char *name) {
	char     cmd[1024 + FILENAME_MAX];
	memset(cmd, 0, 1024 + FILENAME_MAX);

	LayRecordStore   *lof = new LayRecordStore;
	initializeLOF(lof);
	lof->isOutput = false;

	if ((name == NULL) || (strcmp(name, "-") == 0))
		name = NULL;

	errno = 0;
	if        (name == NULL) {
		lof->file = stdin;
	/*
	 * disable zip file support because it does not seem to work on all filesystem with threads
	 } else if (strcasecmp(name+strlen(name)-3, ".gz") == 0) {
		sprintf(cmd, "gzip -dc %s", name);
		lof->file = popen(cmd, "r");
		lof->isPopened = TRUE;
	} else if (strcasecmp(name+strlen(name)-4, ".bz2") == 0) {
		sprintf(cmd, "bzip2 -dc %s", name);
		lof->file = popen(cmd, "r");
		lof->isPopened = TRUE;
	*/} else {
		lof->file       = fopen(name, "r");
	}
	if (errno) {
		fprintf(stderr, "openLayFile()-- Failed to open '%s' for reading: %s\n",
				name, strerror(errno));
		lof = NULL;
	}

	return(lof);
}


void
closeLayFile(LayRecordStore *lof) {

	if (lof == NULL)
		return;

	if (lof->isOutput)
		flushLOF(lof);

	if (lof->isPopened)
		pclose(lof->file);
	else if (lof->file != stdout)
		fclose(lof->file);

	delete[] lof->buffer;
	delete lof;
}

static bool readLayRecord(LayRecordStore *in, LayRecord &r, PBRThreadGlobals *waGlobal, MultiAlignT *ma) {
	clearLayout(r);

	// read in a record
	uint32 count = 0;

	if (in->bufferPos + LAY_HEADER_SIZE > in->bufferLen) {
		uint32 leftovers = in->bufferLen - in->bufferPos;
		for (int i = 0; i < leftovers; i++) {
			in->buffer[i] = in->buffer[in->bufferPos+i];
		}
		in->bufferPos = 0;
		in->bufferLen = AS_UTL_safeRead(in->file,
				in->buffer + leftovers,
				"AS_PBR_readLayRecord",
				sizeof(int32),
				in->bufferMax - leftovers) + leftovers;
	}

	if (in->bufferPos + LAY_HEADER_SIZE > in->bufferLen)
		return(FALSE);

	r.iid = in->buffer[in->bufferPos++];
	count = in->buffer[in->bufferPos++];

	if (ma != NULL) {
		ma->maID = r.iid;
		ma->data.num_frags = count;
		ResetToRange_VA(ma->f_list, ma->data.num_frags);
	}
	for (uint32 iter = 0; iter < count; iter++) {
		OverlapPos o;
		SeqInterval bclr;
		OVSoverlap olap;

		if (in->bufferPos + LAY_RECORD_SIZE > in->bufferLen) {
			uint32 leftovers = in->bufferLen - in->bufferPos;
			for (int i = 0; i < leftovers; i++) {
				in->buffer[i] = in->buffer[in->bufferPos+i];
			}
			in->bufferPos = 0;
			in->bufferLen = AS_UTL_safeRead(in->file,
					in->buffer + leftovers,
					"AS_PBR_readLayRecord",
					sizeof(int32),
					in->bufferMax - leftovers) + leftovers;
		}
		o.ident = in->buffer[in->bufferPos++];
		o.position.bgn = in->buffer[in->bufferPos++];
		o.position.end = in->buffer[in->bufferPos++];
		bclr.bgn = in->buffer[in->bufferPos++];
		bclr.end = in->buffer[in->bufferPos++];
		memcpy(olap.dat.dat, in->buffer + in->bufferPos, sizeof(int32)*AS_OVS_NWORDS);
		in->bufferPos += AS_OVS_NWORDS;

		olap.a_iid = r.iid;
		olap.b_iid = o.ident;
		r.mp.push_back(o);
		r.bClrs[o.ident] = bclr;
		r.bOvls[o.ident] = olap;
		if (r.length < MAX(o.position.bgn, o.position.end)) {
			r.length = MAX(o.position.bgn, o.position.end);
		}

		if (waGlobal != NULL) {
			AS_IID libID = waGlobal->frgToLib[o.ident];
			if (waGlobal->libToOrientation[libID] != AS_READ_ORIENT_UNKNOWN) {
				// save the beginning position of each sequence
				// we can save the beginning position because if the read is fwd it looks like:
				//	PacRead: -------------------------->
				//	fwd read:		--->
				//	rev read:				<---
				// so the begin is always the outer-most positions of the sequences
				r.readStarts[o.ident] = o.position.bgn;
				r.readOri[o.ident] = (o.position.bgn < o.position.end);
			}
		}
		if (ma != NULL) {
			IntMultiPos  *imp = GetIntMultiPos(ma->f_list, iter);
			imp->ident        = o.ident;
			imp->contained    = false;
			imp->parent       = 0;
			imp->ahang        = 0;
			imp->bhang        = 0;
			imp->position.bgn = o.position.bgn;
			imp->position.end = o.position.end;
		}
	}

	return true;
}

bool readLayRecord(LayRecordStore *in, LayRecord &r) {
	return readLayRecord(in, r, NULL, NULL);
}

bool readLayRecord(LayRecordStore *in, LayRecord &r, PBRThreadGlobals *waGlobal) {
	return readLayRecord(in, r, waGlobal, NULL);
}

bool readLayRecord(LayRecordStore *in, LayRecord &r, MultiAlignT *ma) {
	return readLayRecord(in, r, NULL, ma);
}

uint32 writeLayRecord(LayRecordStore *out, LayRecord &r, boost::dynamic_bitset<> *bits, double percentageToStoreInBits, map<AS_IID, uint8> &counts, map<AS_IID, uint32> &largeCounts, double threshold, map<AS_IID, bool> *subset, BinaryOverlapFile *bof) {
	bool isRecordingGaps = false;
	assert(out->isOutput == TRUE);

	if (out->bufferLen + LAY_HEADER_SIZE >= out->bufferMax) {
		AS_UTL_safeWrite(out->file, out->buffer, "AS_PBR_writeLayRecord", sizeof(int32), out->bufferLen);
		out->bufferLen = 0;
	}

	uint32 count = r.mp.size();
	uint32 gapBP = 0;
	uint32 lastEnd = 0;
	uint32 lastGapEnd = 0;
	uint32 lastRecorded = 0;
	uint32 lastRecordedIndex = 0;
	uint32 lastStored = 0;

	if (subset != NULL) {
		count = subset->size();
	}
	out->buffer[out->bufferLen++] = r.iid;
	out->buffer[out->bufferLen++] = count;

	if (VERBOSE) fprintf(stderr, "******Storing long read %d\n", r.iid);
	for (vector<OverlapPos>::const_iterator iter = r.mp.begin(); iter != r.mp.end(); iter++) {
		AS_IID ident = iter->ident;

		if (subset == NULL || (subset != NULL && subset->find(ident) != subset->end())) {
			if (out->bufferLen + LAY_RECORD_SIZE >= out->bufferMax) {
				AS_UTL_safeWrite(out->file, out->buffer, "AS_PBR_writeLayRecord", sizeof(int32), out->bufferLen);
				out->bufferLen = 0;
			}

			// when we have gaps in the tilling, record those sequences that surround the gaps
			// first we keep a running tally of the requested length of surrounding sequences (randomly sampled) at all times in case we see a gap
            if (lastRecordedIndex == 0 && lastRecorded == 0) {
                lastRecordedIndex = iter-r.mp.begin();
                lastRecorded = MIN(iter->position.bgn, iter->position.end);
            }
			if (MIN(iter->position.bgn, iter->position.end) - lastRecorded > MIN_DIST_TO_RECRUIT) {
				// pop one and add new one
			    if (VERBOSE) fprintf(stderr, "Tracking on position %td and unset %d (read %d). Last position I stored is %d, current pos is %d ", iter-r.mp.begin(), lastRecordedIndex, r.mp[lastRecordedIndex].ident, lastRecorded, MIN(iter->position.bgn, iter->position.end));

				bits->set(r.mp[lastRecordedIndex].ident, false);
				do {
				    lastRecordedIndex++;
				} while (subset != NULL && subset->find(r.mp[lastRecordedIndex].ident) == subset->end());
				lastRecorded = MIN(r.mp[lastRecordedIndex].position.bgn, r.mp[lastRecordedIndex].position.end);
				if (VERBOSE) fprintf(stderr, " and shifted to position %d offset %d\n", lastRecordedIndex, lastRecorded);
			}
			bits->set(iter->ident, (drand48() <= percentageToStoreInBits));

			// when we do see a gap, we record the surround sequences of that gap and reset our preceeding range information to start a new range
			// also when we see a sequence that looks repetitive, record its surroundings
            uint32 count = counts[iter->ident];
			if (count == MAX_COV) {
			    count = largeCounts[iter->ident];
            }
			if ((lastEnd != 0 && lastEnd <= MIN(iter->position.bgn, iter->position.end)) ||
			    count > threshold) {
				gapBP += (count > threshold ? MAX(iter->position.bgn, iter->position.end) - MIN(iter->position.bgn, iter->position.end) : MIN(iter->position.bgn, iter->position.end) - lastEnd);
				lastGapEnd = MIN(iter->position.bgn, iter->position.end);
				lastRecordedIndex = lastRecorded = 0;
				isRecordingGaps = true;

				// always record the surrounding reads of a gap
				bits->set(r.mp[lastStored].ident, true);
				bits->set(iter->ident, true);
				if (VERBOSE) fprintf(stderr, "Stored surrounding reads of gap between index %td and surrounding reads: %d and %d in %d\n", iter-r.mp.begin(), r.mp[lastStored].ident, iter->ident, r.iid);
			}

			// we also record a subset of the sequences following a gap
			if (isRecordingGaps) {
                if (VERBOSE) fprintf(stderr, "Setting after on position %td of "F_SIZE_T" (read %d) and set position %d due to gap at %d\n", iter-r.mp.begin(), r.mp.size(), iter->ident, (MAX(iter->position.bgn, iter->position.end)), lastGapEnd);

                if (bits->test(iter->ident) == false) {
                    // drop with random chance
                    bits->set(iter->ident, (drand48() <= percentageToStoreInBits));
                }
                isRecordingGaps = (MAX(iter->position.bgn, iter->position.end) - lastGapEnd) < MIN_DIST_TO_RECRUIT;
			}

			out->buffer[out->bufferLen++] = iter->ident;
			out->buffer[out->bufferLen++] = iter->position.bgn;
			out->buffer[out->bufferLen++] = iter->position.end;
			out->buffer[out->bufferLen++] = r.bClrs[iter->ident].bgn;
			out->buffer[out->bufferLen++] = r.bClrs[iter->ident].end;
			OVSoverlap tmp = r.bOvls[iter->ident];
			SeqInterval pos = iter->position;
			memcpy(out->buffer + out->bufferLen, r.bOvls[iter->ident].dat.dat, sizeof(int32)*AS_OVS_NWORDS);
			out->bufferLen += AS_OVS_NWORDS;

			if (bof != NULL) {
				AS_OVS_writeOverlap(bof, &r.bOvls[iter->ident]);
			}

			if (lastEnd < MAX(iter->position.bgn, iter->position.end)) {
				lastEnd = MAX(iter->position.bgn, iter->position.end);
			}
            lastStored = iter-r.mp.begin();
		}
	}

	if (gapBP == 0) {
		bits->reset();
	} else if (VERBOSE) {
	    fprintf(stderr, "******For long read %d the following are set\n", r.iid);
	    for (int i = 0; i < bits->size(); i++) {
	        if (bits->test(i)) {
	            fprintf(stderr, "Short read ID %d\n", i);
	        }
	    }
	    fprintf(stderr, "************************************************\n");
	}
	return gapBP;
}
