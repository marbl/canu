
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-NOV-26 to 2015-AUG-10
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "gkStore.H"

#include "AS_UTL_fileIO.H"
#include "AS_UTL_alloc.H"


gkStore *gkStore::_instance      = NULL;
uint32   gkStore::_instanceCount = 0;



bool
gkRead::gkRead_loadData(gkReadData *readData, void *blobs) {

  readData->_read = this;

  //  The resize will only increase the space.  if the new is less than the max, it returns immediately.

  resizeArrayPair(readData->_seq, readData->_qlt, readData->_seqAlloc, readData->_seqAlloc, (uint32)_seqLen+1, resizeArray_doNothing);

  //  Where, or where!, is the data?

  uint64  offset = _mPtr;

  //  One might be tempted to set the readData blob to point to the blob data in the mmap,
  //  but doing so will cause it to be written out again.

  readData->_blobLen = 0;
  readData->_blobMax = 0;
  readData->_blob    = NULL;

  //  Instead, we'll use someting horribly similar.

  uint8  *blob    = ((uint8 *)blobs) + offset;
  char    chunk[5];

  if ((blob[0] != 'B') && (blob[1] != 'L') && (blob[2] != 'O') && (blob[3] != 'B'))
    fprintf(stderr, "Index error in read "F_U32" %c mPtr "F_U64" pID "F_U64" expected BLOB, got %c%c%c%c\n",
            gkRead_readID(),
            '?', //(_numberOfPartitions == 0) ? 'm' : 'p',
            _mPtr, _pID,
            blob[0], blob[1], blob[2], blob[3]);
  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  uint32  blobLen = *((uint32 *)blob + 1);

  //fprintf(stderr, "BLOB len %u\n", blobLen);

  blob += 8;

  while ((blob[0] != 'S') ||
         (blob[1] != 'T') ||
         (blob[2] != 'O') ||
         (blob[3] != 'P')) {
    chunk[0] = blob[0];
    chunk[1] = blob[1];
    chunk[2] = blob[2];
    chunk[3] = blob[3];
    chunk[4] = 0;

    uint32   chunkLen = *((uint32 *)blob + 1);

#if 1
    // Fix for stores built between 9/3/15-9/7/15 when QVs were changed to uniform value but nothing was stored on disk
    //    Loading from these stores lead to random uninitialized QVs
    //    set all QVs to a default value
    // Should be unnecessary after Dec 8, 2015 and can be removed
    for (uint32 ii=0; ii<_seqLen; ii++)
      readData->_qlt[ii] = 20;
#endif

    //fprintf(stderr, "%s len %u\n", chunk, chunkLen);

    if      (strncmp(chunk, "VERS", 4) == 0) {
    }

    else if (strncmp(chunk, "QSEQ", 4) == 0) {
      //fprintf(stderr, "QSEQ not supported.\n");
    }

    else if (strncmp(chunk, "USEQ", 4) == 0) {
      assert(_seqLen <= chunkLen);
      assert(_seqLen <= readData->_seqAlloc);
      memcpy(readData->_seq, blob + 8, _seqLen);
      readData->_seq[_seqLen] = 0;
    }

    else if (strncmp(chunk, "UQLT", 4) == 0) {
      assert(_seqLen <= chunkLen);
      assert(_seqLen <= readData->_seqAlloc);
      memcpy(readData->_qlt, blob + 8, _seqLen);
      readData->_qlt[_seqLen] = 0;

#if 1
      //  Fix for older gkpStore that encoded QV's with offset '0'.  Wasn't RIFF format supposed to solve problems like this?
      //    Old encoding is ASCII was from '0' = 48 to 'l' = 108.
      //    New encoding is integer from 0 to 60.
      //  So, if we see:
      //    0-47,  we're new format.
      //    48-60  we're either (but STRONGLY likely to be old).
      //    61-108 we're old format.
      //  Most reads are well below qv=40, so this will be easy.

      bool  isOld = false;

      for (uint32 ii=0; ii<_seqLen; ii++) {
        if (readData->_qlt[ii] < 48) {
          isOld = false;
          break;
        }

        else if (readData->_qlt[ii] < 61) {
          isOld = true;
        }

        else {
          isOld = true;
          break;
        }
      }

      if (isOld)
        for (uint32 ii=0; ii<_seqLen; ii++)
          readData->_qlt[ii] -= '0';
#endif
    }

    else if (strncmp(chunk, "2SEQ", 4) == 0) {
      gkRead_decode2bit(blob + 8, chunkLen, readData->_seq, _seqLen);
    }

    else if (strncmp(chunk, "3SEQ", 4) == 0) {
      gkRead_decode3bit(blob + 8, chunkLen, readData->_seq, _seqLen);
    }

    else if (strncmp(chunk, "4QLT", 4) == 0) {
      gkRead_decode4bit(blob + 8, chunkLen, readData->_qlt, _seqLen);
    }

    else if (strncmp(chunk, "5QLT", 4) == 0) {
      gkRead_decode5bit(blob + 8, chunkLen, readData->_qlt, _seqLen);
    }

    else if (strncmp(chunk, "QVAL", 4) == 0) {
      uint32  qval = *((uint32 *)blob + 2);

      for (uint32 ii=0; ii<_seqLen; ii++)
        readData->_qlt[ii] = qval;
    }

    else {
      fprintf(stderr, "gkRead::gkRead_loadData()--  unknown chunk type '%s' skipped\n", chunk);
    }

    blob += 4 + 4 + chunkLen;
  }

  return(true);
};




//  Dump a block of encoded data to disk, then update the gkRead to point to it.
//
void
gkStore::gkStore_stashReadData(gkRead *read, gkReadData *data) {

  assert(_blobsFile != NULL);

  read->_mPtr = AS_UTL_ftell(_blobsFile);
  read->_pID  = _partitionID;                //  0 if not partitioned

  //fprintf(stderr, "STASH read %u at position "F_SIZE_T"\n", read->gkRead_readID(), AS_UTL_ftell(_blobsFile));

  AS_UTL_safeWrite(_blobsFile,
                   data->_blob,
                   "gkStore_stashReadData::blob",
                   sizeof(char),
                   data->_blobLen);
}



//  Load read metadata and data from a stream.
//
void
gkStore::gkStore_loadReadFromStream(FILE *S, gkRead *read, gkReadData *readData) {
  char    tag[5];
  uint32  size;

  //  Mark this as a read.  Needed for tgTig::loadFromStreamOrLayout(), and loading this stuff in
  //  utgcns.

  AS_UTL_safeRead(S, tag, "gkStore::gkStore_loadReadFromStream::tag", sizeof(char), 4);

  if (strncmp(tag, "READ", 4) != 0)
    fprintf(stderr, "Failed to load gkRead, got tag '%c%c%c%c' (0x%02x 0x%02x 0x%02x 0x%02x), expected 'READ'.\n",
            tag[0], tag[1], tag[2], tag[3],
            tag[0], tag[1], tag[2], tag[3]), exit(1);

  //  Load the read metadata

  AS_UTL_safeRead(S, read, "gkStore::gkStore_loadReadFromStream::read", sizeof(gkRead), 1);

  //  With some pain, we read the BLOB and its length, then allocate space for the blob
  //  and finsh reading it.

  AS_UTL_safeRead(S,  tag,  "gkStore::gkStore_loadReadFromStream::blob", sizeof(int8),   4);
  AS_UTL_safeRead(S, &size, "gkStore::gkStore_loadReadFromStream::size", sizeof(uint32), 1);

  uint8 *blob = new uint8 [8 + size];

  memcpy(blob,    tag,  sizeof(uint8)  * 4);
  memcpy(blob+4, &size, sizeof(uint32) * 1);

  AS_UTL_safeRead(S, blob+8, "gkStore::gkStore_loadReadFromStream::blob", sizeof(char), size);

  //  Unpack the blob into a readData

  read->_mPtr = 0;
  read->gkRead_loadData(readData, blob);

  //  And, that's it!  Sweet!
}


//  Dump the read metadata and read data to a stream.
//
void
gkStore::gkStore_saveReadToStream(FILE *S, uint32 id) {

  //  Mark this as a read.  Needed for tgTig::loadFromStreamOrLayout(), and loading this stuff in
  //  utgcns.

  fprintf(S, "READ");

  //  Dump the read metadata

  gkRead  *read = gkStore_getRead(id);

  AS_UTL_safeWrite(S, read, "gkStore::gkStore_saveReadToStream::read", sizeof(gkRead), 1);

  //  Figure out where the blob actually is, and make sure that it really is a blob

  uint8  *blob    = (uint8 *)_blobs + read->_mPtr;
  uint32  blobLen = 8 + *((uint32 *)blob + 1);

  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  //  Write the blob to the stream

  AS_UTL_safeWrite(S, blob, "gkStore::gkStore_saveReadToStream::blob", sizeof(char), blobLen);
}




//  Store the 'len' bytes of data in 'dat' into the class-managed _blob data block.
//  Ensures that the _blob block is appropriately padded to maintain 32-bit alignment.
//
void
gkReadData::gkReadData_encodeBlobChunk(char const *tag,
                                       uint32      len,
                                       void       *dat) {

  //  Allocate an initial blob if we don't have one

  if (_blobMax == 0) {
    _blobLen = 0;
    _blobMax = 1048576;
    _blob    = new uint8 [_blobMax];
  }

  //  Or make it bigger

  while (_blobMax + 8 + len < _blobMax) {
    _blobMax *= 2;
    uint8 *b  = new uint8 [_blobMax];
    memcpy(b, _blob, sizeof(uint8) * _blobLen);
    delete [] _blob;
    _blob = b;
  }

  //  Figure out how much padding we need to add

  uint32 pad = 4 - (len % 4);

  if (pad == 4)
    pad = 0;

  //  Copy in the chunk id and padded length

  len += pad;

  memcpy(_blob + _blobLen,  tag, sizeof(char) * 4);     _blobLen += sizeof(char) * 4;
  memcpy(_blob + _blobLen, &len, sizeof(uint32));       _blobLen += sizeof(uint32);

  len -= pad;

  //  Then the unpadded data and any padding.

  memcpy(_blob + _blobLen,  dat, sizeof(uint8) * len);  _blobLen += sizeof(uint8) * len;

  if (pad > 2)  _blob[_blobLen++] = 0;
  if (pad > 1)  _blob[_blobLen++] = 0;
  if (pad > 0)  _blob[_blobLen++] = 0;

  //  Finally, update the total blob length.

  _blobLen -= 8;

  memcpy(_blob + 4, &_blobLen, sizeof(uint32));

  _blobLen += 8;
}



gkReadData *
gkRead::gkRead_encodeSeqQlt(char *H, char *S, char *Q, uint32 qv) {
  gkReadData *rd = new gkReadData;

  uint32  RID  = _readID;    //  Debugging

  //  If there is a QV string, ensure that the lengths are the same.  If not, trim or pad the QVs.
  //  Then, convert the expected Sanger-encoded QV's (base='!') to be just integers.

  uint32  Slen = _seqLen = strlen(S);
  uint32  Qlen = 0;

  if (Q[0] != 0) {
    Qlen = strlen(Q);

    if (Slen < Qlen) {
      fprintf(stderr, "-- WARNING:  read '%s' sequence length %u != quality length %u; quality bases truncated.\n",
              H, Slen, Qlen);
      Q[Slen] = 0;
    }

    if (Slen > Qlen) {
      fprintf(stderr, "-- WARNING:  read '%s' sequence length %u != quality length %u; quality bases padded.\n",
              H, Slen, Qlen);
      for (uint32 ii=Qlen; ii<Slen; ii++)
        Q[ii] = Q[Qlen-1];
    }

    for (uint32 ii=0; ii<Qlen; ii++)
      Q[ii] -= '!';
  }

  //  Compute the preferred encodings.  If either fail, the length is set to zero, and ...

  uint8   *seq = NULL;
  uint8   *qlt = NULL;

  uint32  seq2Len = gkRead_encode2bit(seq, S, Slen);
  uint32  seq3Len = 0;
  uint32  qlt4Len = gkRead_encode4bit(qlt, Q, Qlen);
  uint32  qlt5Len = 0;

  //  ... the non-preferred encoding will be computed.  If this too fails, sequences/qualities will
  //  be stored unencoded.

  if (seq2Len == 0)
    seq3Len = gkRead_encode3bit(seq, S, Slen);

  if (qlt4Len == 0)
    qlt5Len = gkRead_encode5bit(qlt, Q, Qlen);

  //  Encode the data into chunks in the blob.

  uint32  blobVers = 0x00000001;

  rd->gkReadData_encodeBlobChunk("BLOB",       0,  NULL);
  rd->gkReadData_encodeBlobChunk("VERS",       4, &blobVers);

  if      (seq2Len > 0)
    rd->gkReadData_encodeBlobChunk("2SEQ", seq2Len, seq);    //  Two-bit encoded sequence (ACGT only)
  else if (seq3Len > 0)
    rd->gkReadData_encodeBlobChunk("3SEQ", seq3Len, seq);    //  Three-bit encoded sequence (ACGTN)
  else
    rd->gkReadData_encodeBlobChunk("USEQ", Slen, S);         //  Unencoded sequence

  if      (qlt4Len > 0)
    rd->gkReadData_encodeBlobChunk("4QLT", qlt4Len, qlt);    //  Four-bit (0-15) encoded QVs
  else if (qlt5Len > 0)
    rd->gkReadData_encodeBlobChunk("5QLT", qlt5Len, qlt);    //  Five-bit (0-32) encoded QVs
  else if (Q[0] == 0)
    rd->gkReadData_encodeBlobChunk("QVAL", 4, &qv);          //  Constant QV for every base
  else
    rd->gkReadData_encodeBlobChunk("UQLT", Qlen, Q);         //  Unencoded quality

  rd->gkReadData_encodeBlobChunk("STOP", 0,  NULL);

  //  Cleanup.  Restore the QV's.  Delete temporary storage.

  for (uint32 ii=0; ii<Qlen; ii++)
    Q[ii] += '!';

  delete [] seq;
  delete [] qlt;

  return(rd);
}



#if 0
//  Not implemented.
gkReadData *
gkRead::gkRead_encodePacBio(char *H, char *S, char *Q) {
  gkReadData *rd = new gkReadData;

  return(rd);
}

//  Not implemented.
gkReadData *
gkRead::gkRead_encodeMinION(char *H, char *S, char *Q) {
  gkReadData *rd = new gkReadData;

  return(rd);
}
#endif



////////////////////////////////////////
//
//  gkLibrary is lightweight, except for three functions that need to parse strings
//

void
gkLibrary::gkLibrary_parsePreset(char *p) {

  if (strcasecmp(p, "contig") == 0) {
    _readCorrection   = GK_CORRECTION_NONE;
    _readType         = GK_READTYPE_CONTIG;
    return;
  }

  if (strcasecmp(p, "pacbio-raw") == 0) {
    _readCorrection   = GK_CORRECTION_CONSENSUS;
    _readType         = GK_READTYPE_PACBIO_RAW;
    _checkForSubReads = true;
    return;
  }

  if (strcasecmp(p, "pacbio-corrected") == 0) {
    _readCorrection   = GK_CORRECTION_NONE;
    _readType         = GK_READTYPE_PACBIO_CORRECTED;
    return;
  }

  if (strcasecmp(p, "nanopore-raw") == 0) {
    _readCorrection   = GK_CORRECTION_CONSENSUS;
    _readType         = GK_READTYPE_NANOPORE_RAW;
    _checkForSubReads = true;
    return;
  }

  if (strcasecmp(p, "nanopore-corrected") == 0) {
    _readCorrection   = GK_CORRECTION_NONE;
    _readType         = GK_READTYPE_NANOPORE_CORRECTED;
    return;
  }

  fprintf(stderr, "gkLibrary::gkLibrary_parsePreset()--  ERROR: unknown preset '%s'\n", p);
  exit(1);
}


void
gkLibrary::gkLibrary_setReadType(char *t) {

  if      (strcasecmp(t, "generic") == 0)
    _readType = GK_READTYPE_GENERIC;

  else if (strcasecmp(t, "contig") == 0)
    _readType = GK_READTYPE_CONTIG;

  else if (strcasecmp(t, "pacbio_raw") == 0)
    _readType = GK_READTYPE_PACBIO_RAW;

  else if (strcasecmp(t, "pacbio_corrected") == 0)
    _readType = GK_READTYPE_PACBIO_CORRECTED;

  else if (strcasecmp(t, "nanopore_raw") == 0)
    _readType = GK_READTYPE_NANOPORE_RAW;

  else if (strcasecmp(t, "nanopore_corrected") == 0)
    _readType = GK_READTYPE_NANOPORE_CORRECTED;

  else
    fprintf(stderr, "gkLibrary::gkLibrary_setReadType()--  ERROR: unknown read type '%s'\n",
            t), exit(1);
}


char const *
gkLibrary::gkLibrary_readTypeString(void) {
  switch (_readType) {
    case GK_READTYPE_GENERIC:             return("generic");             break;
    case GK_READTYPE_CONTIG:              return("contig");              break;
    case GK_READTYPE_PACBIO_RAW:          return("pacbio-raw");          break;
    case GK_READTYPE_PACBIO_CORRECTED:    return("pacbio-corrected");    break;
    case GK_READTYPE_NANOPORE_RAW:        return("nanopore-raw");        break;
    case GK_READTYPE_NANOPORE_CORRECTED:  return("nanopore-corrected");  break;
    default:                              return("invalid");             break;
  }
}


void
gkLibrary::gkLibrary_setReadCorrection(char *t) {

  if      (strcasecmp(t, "none") == 0)
    _finalTrim = GK_CORRECTION_NONE;

  else if (strcasecmp(t, "consensus") == 0)
    _finalTrim = GK_CORRECTION_CONSENSUS;

  else if (strcasecmp(t, "mer") == 0)
    _finalTrim = GK_CORRECTION_MER;

  else
    fprintf(stderr, "gkLibrary::gkLibrary_setReadCorrection()--  ERROR: unknown read correction '%s'\n",
            t), exit(1);
}


char const *
gkLibrary::gkLibrary_readCorrectionString(void) {
  switch (_readCorrection) {
    case GK_CORRECTION_NONE:              return("none");                break;
    case GK_CORRECTION_CONSENSUS:         return("consensus");           break;
    case GK_CORRECTION_MER:               return("mer");                 break;
    default:                              return("invalid");             break;
  }
}


void
gkLibrary::gkLibrary_setFinalTrim(char *t) {

  if      (strcasecmp(t, "none") == 0)
    _finalTrim = GK_FINALTRIM_NONE;

  else if (strcasecmp(t, "largest") == 0)
    _finalTrim = GK_FINALTRIM_LARGEST_COVERED;

  else if (strcasecmp(t, "bestedge") == 0)
    _finalTrim = GK_FINALTRIM_BEST_EDGE;

  else
    fprintf(stderr, "gkLibrary::gkLibrary_setFinalTrim()--  ERROR: unknown final trim '%s'\n",
            t), exit(1);
}


char const *
gkLibrary::gkLibrary_finalTrimString(void) {
  switch (_finalTrim) {
    case GK_FINALTRIM_NONE:               return("none");                break;
    case GK_FINALTRIM_LARGEST_COVERED:    return("largestCovered");      break;
    case GK_FINALTRIM_BEST_EDGE:          return("bestEdge");            break;
    default:                              return("invalid");             break;
  }
}








////////////////////////////////////////

//  The N valid modes for a 'new gkpStore' call:
//
//  1)  Add new reads/libraries, modify old ones.  gkStore(path, true, true)
//  2)  No addition, but can modify old ones.      gkStore(path, true)
//  3)  No addition, no modification.              gkStore(path);
//
gkStore::gkStore(char const *path, gkStore_mode mode, uint32 partID) {
  char    name[FILENAME_MAX + 5];

  memset(_storePath, 0, sizeof(char) * FILENAME_MAX);
  memset(_storeName, 0, sizeof(char) * FILENAME_MAX);

  strncpy(_storePath, path, FILENAME_MAX-1);
  strncpy(_storeName, path, FILENAME_MAX-1);  //  Broken.

  sprintf(name, "%s/info", _storePath);

  //  If the info file exists, load it.

  if (AS_UTL_fileExists(name, false, false) == true) {
    errno = 0;
    FILE *I = fopen(name, "r");
    AS_UTL_safeRead(I, &_info, "gkStore::_info", sizeof(gkStoreInfo), 1);
    fclose(I);
  }

  //  Check sizes are correct.

  uint32  failed = 0;

  if (_info.gkLibrarySize      != sizeof(gkLibrary))
    failed += fprintf(stderr, "ERROR:  gkLibrary size in store = "F_U32", differs from executable = "F_SIZE_T"\n",
                      _info.gkLibrarySize, sizeof(gkLibrary));

  if (_info.gkReadSize         != sizeof(gkRead))
    failed += fprintf(stderr, "ERROR:  gkRead size in store = "F_U32", differs from executable = "F_SIZE_T"\n",
                      _info.gkReadSize, sizeof(gkRead));

  if (_info.gkMaxLibrariesBits != AS_MAX_LIBRARIES_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_LIBRARIES_BITS in store = "F_U32", differs from executable = "F_U32"\n",
                      _info.gkMaxLibrariesBits, AS_MAX_LIBRARIES_BITS);

  if (_info.gkLibraryNameSize  != LIBRARY_NAME_SIZE)
    failed += fprintf(stderr, "ERROR:  LIBRARY_NAME_SIZE in store = "F_U32", differs from executable = "F_U32"\n",
                      _info.gkLibraryNameSize, LIBRARY_NAME_SIZE);

  if (_info.gkMaxReadBits      != AS_MAX_READS_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_READS_BITS in store = "F_U32", differs from executable = "F_U32"\n",
                      _info.gkMaxReadBits, AS_MAX_READS_BITS);

  if (_info.gkMaxReadLenBits   != AS_MAX_READLEN_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_READLEN_BITS in store = "F_U32", differs from executable = "F_U32"\n",
                      _info.gkMaxReadLenBits, AS_MAX_READLEN_BITS);

  if (failed)
    fprintf(stderr, "ERROR:\nERROR:  Can't open store '%s': parameters in src/AS_global.H are incompatible with the store.\n", _storePath), exit(1);

  assert(_info.gkLibrarySize      == sizeof(gkLibrary));
  assert(_info.gkReadSize         == sizeof(gkRead));

  assert(_info.gkMaxLibrariesBits == AS_MAX_LIBRARIES_BITS);
  assert(_info.gkLibraryNameSize  == LIBRARY_NAME_SIZE);
  assert(_info.gkMaxReadBits      == AS_MAX_READS_BITS);
  assert(_info.gkMaxReadLenBits   == AS_MAX_READLEN_BITS);

  //  Clear ourself, to make valgrind happier.

  _librariesMMap          = NULL;
  _librariesAlloc         = 0;
  _libraries              = NULL;

  _readsMMap              = NULL;
  _readsAlloc             = 0;
  _reads                  = NULL;

  _blobsMMap              = NULL;
  _blobs                  = NULL;
  _blobsFile              = NULL;

  _mode                   = mode;

  _numberOfPartitions     = 0;
  _partitionID            = 0;
  _readIDtoPartitionIdx   = NULL;
  _readIDtoPartitionID    = NULL;
  _readsPerPartition      = NULL;
  //_readsInThisPartition   = NULL;

  //
  //  READ ONLY
  //

  if ((mode == gkStore_readOnly) &&
      (partID == UINT32_MAX)) {
    //fprintf(stderr, "gkStore()--  opening '%s' for read-only access.\n", _storePath);

    if (AS_UTL_fileExists(_storePath, true, false) == false) {
      fprintf(stderr, "gkStore()--  failed to open '%s' for read-only access: store doesn't exist.\n", _storePath);
      exit(1);
    }

    sprintf(name, "%s/libraries", _storePath);
    _librariesMMap = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _libraries     = (gkLibrary *)_librariesMMap->get(0);

    sprintf(name, "%s/reads", _storePath);
    _readsMMap     = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _reads         = (gkRead *)_readsMMap->get(0);

    sprintf(name, "%s/blobs", _storePath);
    _blobsMMap     = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _blobs         = (void *)_blobsMMap->get(0);
  }

  //
  //  MODIFY, NO APPEND (also for building a partitioned store)
  //

  else if ((mode == gkStore_modify) &&
           (partID == UINT32_MAX)) {
    //fprintf(stderr, "gkStore()--  opening '%s' for read-write access.\n", _storePath);

    if (AS_UTL_fileExists(_storePath, true, false) == false) {
      fprintf(stderr, "gkStore()--  failed to open '%s' for read-write access: store doesn't exist.\n", _storePath);
      exit(1);
    }

    sprintf(name, "%s/libraries", _storePath);
    _librariesMMap = new memoryMappedFile (name, memoryMappedFile_readWrite);
    _libraries     = (gkLibrary *)_librariesMMap->get(0);

    sprintf(name, "%s/reads", _storePath);
    _readsMMap     = new memoryMappedFile (name, memoryMappedFile_readWrite);
    _reads         = (gkRead *)_readsMMap->get(0);

    sprintf(name, "%s/blobs", _storePath);
    _blobsMMap     = new memoryMappedFile (name, memoryMappedFile_readWrite);
    _blobs         = (void *)_blobsMMap->get(0);
  }

  //
  //  MODIFY, APPEND, open mmap'd files, but copy them entirely to local memory
  //

  else if ((mode == gkStore_extend) &&
           (partID == UINT32_MAX)) {
    //fprintf(stderr, "gkStore()--  opening '%s' for read-write and append access.\n", _storePath);

    if (AS_UTL_fileExists(_storePath, true, true) == false)
      AS_UTL_mkdir(_storePath);

    _librariesAlloc = MAX(64, 2 * _info.numLibraries);
    _libraries      = new gkLibrary [_librariesAlloc];

    sprintf(name, "%s/libraries", _storePath);
    if (AS_UTL_fileExists(name, false, false) == true) {
      _librariesMMap  = new memoryMappedFile (name, memoryMappedFile_readOnly);

      memcpy(_libraries, _librariesMMap->get(0), sizeof(gkLibrary) * (_info.numLibraries + 1));

      delete _librariesMMap;
      _librariesMMap = NULL;;
    }

    _readsAlloc     = MAX(128, 2 * _info.numReads);
    _reads          = new gkRead [_readsAlloc];

    sprintf(name, "%s/reads", _storePath);
    if (AS_UTL_fileExists(name, false, false) == true) {
      _readsMMap      = new memoryMappedFile (name, memoryMappedFile_readOnly);

      memcpy(_reads, _readsMMap->get(0), sizeof(gkRead) * (_info.numReads + 1));

      delete _readsMMap;
      _readsMMap = NULL;
    }

    sprintf(name, "%s/blobs", _storePath);

    _blobsMMap     = NULL;
    _blobs         = NULL;

    errno = 0;
    _blobsFile     = fopen(name, "a+");
    if (errno)
      fprintf(stderr, "gkStore()--  Failed to open blobs file '%s' for appending: %s\n",
              name, strerror(errno)), exit(1);
  }

  //
  //  PARTITIONED, no modifications, no appends
  //
  //  BIG QUESTION: do we want to partition the read metadata too, or is it small enough
  //  to load in every job?  For now, we load all the metadata.

  else if ((mode == gkStore_readOnly) &&
           (partID != UINT32_MAX)) {
    //fprintf(stderr, "gkStore()--  opening '%s' partition '%u' for read-only access.\n", _storePath, partID);

    //  For partitioned reads, we need to have a uint32 map of readID to partitionReadID so we can
    //  lookup the metadata in the partitoned _reads data.  This is 4 bytes per read, compared to 24
    //  bytes for the full meta data.  Assuming 100x of 3kb read coverage on human, that's 100
    //  million reads, so 0.400 GB vs 2.4 GB.

    sprintf(name, "%s/partitions/map", _storePath);

    errno = 0;
    FILE *F = fopen(name, "r");
    if (errno)
      fprintf(stderr, "gkStore::gkStore()-- failed to open '%s' for reading: %s\n",
              name, strerror(errno)), exit(1);

    AS_UTL_safeRead(F, &_numberOfPartitions, "gkStore::_numberOfPartitions", sizeof(uint32), 1);

    _partitionID            = partID;
    _readsPerPartition      = new uint32 [_numberOfPartitions   + 1];  //  No zeroth element in any of these
    _readIDtoPartitionID    = new uint32 [gkStore_getNumReads() + 1];
    _readIDtoPartitionIdx   = new uint32 [gkStore_getNumReads() + 1];

    AS_UTL_safeRead(F, _readsPerPartition,    "gkStore::_readsPerPartition",    sizeof(uint32), _numberOfPartitions   + 1);
    AS_UTL_safeRead(F, _readIDtoPartitionID,  "gkStore::_readIDtoPartitionID",  sizeof(uint32), gkStore_getNumReads() + 1);
    AS_UTL_safeRead(F, _readIDtoPartitionIdx, "gkStore::_readIDtoPartitionIdx", sizeof(uint32), gkStore_getNumReads() + 1);

    fclose(F);

    sprintf(name, "%s/libraries", _storePath);
    _librariesMMap = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _libraries     = (gkLibrary *)_librariesMMap->get(0);
    //fprintf(stderr, " -- openend '%s' at "F_X64"\n", name, _libraries);

    sprintf(name, "%s/partitions/reads.%04"F_U32P"", _storePath, partID);
    _readsMMap     = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _reads         = (gkRead *)_readsMMap->get(0);
    //fprintf(stderr, " -- openend '%s' at "F_X64"\n", name, _reads);

    sprintf(name, "%s/partitions/blobs.%04"F_U32P"", _storePath, partID);
    _blobsMMap     = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _blobs         = (void *)_blobsMMap->get(0);
    //fprintf(stderr, " -- openend '%s' at "F_X64"\n", name, _blobs);
  }

  //  Info only, no access to reads or libraries.

  else if (mode == gkStore_infoOnly) {
    //fprintf(stderr, "gkStore()--  opening '%s' for info-only access.\n", _storePath);
  }

  else {
    fprintf(stderr, "gkStore::gkStore()-- invalid mode '%s' with partition ID %u.\n",
            toString(mode), partID);
    assert(0);
  }
}






gkStore::~gkStore() {
  char   N[FILENAME_MAX];
  FILE  *F;

  //  Should check that inf on disk is the same as inf in memory, and update if needed.

  bool   needsInfoUpdate = false;

  //  Write N+1 because we write, but don't count, the [0] element.

  if        (_librariesMMap) {
    delete _librariesMMap;

  } else if (_libraries) {
    sprintf(N, "%s/libraries", gkStore_path());
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "gkStore::~gkStore()-- failed to open '%s' for writing: %s\n",
              N, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, _libraries, "libraries", sizeof(gkLibrary), gkStore_getNumLibraries() + 1);
    fclose(F);

    delete [] _libraries;

    needsInfoUpdate = true;
  }


  if        (_readsMMap) {
    delete _readsMMap;

  } else if (_reads) {
    sprintf(N, "%s/reads", gkStore_path());
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "gkStore::~gkStore()-- failed to open '%s' for writing: %s\n",
              N, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, _reads, "reads", sizeof(gkRead), gkStore_getNumReads() + 1);
    fclose(F);

    delete [] _reads;

    needsInfoUpdate = true;
  }


  if (needsInfoUpdate) {
    sprintf(N, "%s/info", gkStore_path());
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "gkStore::~gkStore()-- failed to open '%s' for writing: %s\n",
              N, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, &_info, "info", sizeof(gkStoreInfo), 1);
    fclose(F);


    sprintf(N, "%s/info.txt", gkStore_path());
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "gkStore::~gkStore()-- failed to open '%s' for writing: %s\n",
              N, strerror(errno)), exit(1);

    _info.writeInfoAsText(F);

    fclose(F);
  }


  if (_blobsMMap)
    delete _blobsMMap;

  if (_blobsFile)
    fclose(_blobsFile);

  delete [] _readIDtoPartitionIdx;
  delete [] _readIDtoPartitionID;
  delete [] _readsPerPartition;
};


gkLibrary *
gkStore::gkStore_addEmptyLibrary(char const *name) {

  assert(_librariesMMap == NULL);
  assert(_info.numLibraries <= _librariesAlloc);

  //  Library[0] doesn't exist, see comments in addEmptyRead below.

  _info.numLibraries++;

  if (_info.numLibraries == _librariesAlloc)
    increaseArray(_libraries, _info.numLibraries, _librariesAlloc, 128);

  //  Bullet proof the library name - so we can make files with this prefix.

  char    libname[LIBRARY_NAME_SIZE];
  uint32  libnamepos = 0;
  bool    modified   = false;
  bool    truncated  = false;

  for (char const *orig=name; *orig; orig++) {
    if        (*orig == '/') {
      libname[libnamepos++] = '_';
      libname[libnamepos++] = '-';
      libname[libnamepos++] = '_';
      modified = true;

    } else if (isspace(*orig) == 0) {
      libname[libnamepos++] = *orig;

    } else {
      libname[libnamepos++] = '_';
      modified = true;
    }

    if (libnamepos >= LIBRARY_NAME_SIZE) {
      libname[LIBRARY_NAME_SIZE-1] = 0;
      truncated = true;
      break;
    }
  }

  libname[libnamepos] = 0;

#if 0
  if (modified || truncated)
    fprintf(stderr, "gkStore_addEmptyLibrary()--  added library '%s' (original name '%s')\n",
            libname, name);
  else
    fprintf(stderr, "gkStore_addEmptyLibrary()--  added library '%s'\n",
            libname);
#endif

  _libraries[_info.numLibraries] = gkLibrary();

  strncpy(_libraries[_info.numLibraries]._libraryName, libname, LIBRARY_NAME_SIZE-1);

  _libraries[_info.numLibraries]._libraryID = _info.numLibraries;

  return(_libraries + _info.numLibraries);
}


gkRead *
gkStore::gkStore_addEmptyRead(gkLibrary *lib) {

  assert(_readsMMap == NULL);
  assert(_info.numReads <= _readsAlloc);
  assert(_mode != gkStore_readOnly);

  //  We reserve the zeroth read for "null".  This is easy to accomplish
  //  here, just pre-increment the number of reads.  However, we need to be sure
  //  to iterate up to and including _info.numReads.

  _info.numReads++;

  if (_info.numReads == _readsAlloc)
    increaseArray(_reads, _info.numReads, _readsAlloc, _info.numReads/2);

  _reads[_info.numReads]            = gkRead();
  _reads[_info.numReads]._readID    = _info.numReads;
  _reads[_info.numReads]._libraryID = lib->gkLibrary_libraryID();

  //fprintf(stderr, "ADD READ %u = %u alloc = %u\n", _info.numReads, _reads[_info.numReads]._readID, _readsAlloc);

  return(_reads + _info.numReads);
}





void
gkRead::gkRead_copyDataToPartition(void *blobs, FILE **partfiles, uint64 *partfileslen, uint32 partID) {

  //  Stash away the location of the partitioned data

  assert(partfileslen[partID] == AS_UTL_ftell(partfiles[partID]));

  //  Figure out where the blob actually is, and make sure that it really is a blob

  uint8  *blob    = (uint8 *)blobs + _mPtr;
  uint32  blobLen = 8 + *((uint32 *)blob + 1);

  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  //  Write the blob to the partition, update the length of the partition

  AS_UTL_safeWrite(partfiles[partID], blob, "gkRead::gkRead_copyDataToPartition::blob", sizeof(char), blobLen);

  //  Update the read to the new location of the blob in the partitioned data.

  _mPtr = partfileslen[partID];
  _pID  = partID;

  //  And finalize by remembering the length.

  partfileslen[partID] += blobLen;

  assert(partfileslen[partID] == AS_UTL_ftell(partfiles[partID]));

}


void
gkStore::gkStore_buildPartitions(uint32 *partitionMap) {
  char              name[FILENAME_MAX];

  //  Store cannot be partitioned already, and it must be readOnly (for safety) as we don't need to
  //  be changing any of the normal store data.

  assert(_numberOfPartitions == 0);
  assert(_mode               == gkStore_readOnly);

  //  Figure out what the last partition is

  uint32  maxPartition = 0;
  uint32  unPartitioned = 0;

  assert(partitionMap[0] == UINT32_MAX);

  for (uint32 fi=1; fi<=gkStore_getNumReads(); fi++) {
    if (partitionMap[fi] == UINT32_MAX)
      unPartitioned++;

    else if (maxPartition < partitionMap[fi])
      maxPartition = partitionMap[fi];
  }

  fprintf(stderr, "Found "F_U32" unpartitioned reads and maximum partition of "F_U32"\n",
          unPartitioned, maxPartition);

  //  Create the partitions by opening N copies of the data stores,
  //  and writing data to each.

  FILE         **blobfiles    = new FILE * [maxPartition + 1];
  uint64        *blobfileslen = new uint64 [maxPartition + 1];            //  Offset, in bytes, into the blobs file
  FILE         **readfiles    = new FILE * [maxPartition + 1];
  uint32        *readfileslen = new uint32 [maxPartition + 1];            //  aka _readsPerPartition
  uint32        *readIDmap    = new uint32 [gkStore_getNumReads() + 1];   //  aka _readIDtoPartitionIdx

  //  Be nice and put all the partitions in a subdirectory.

  sprintf(name,"%s/partitions", _storePath);

  if (AS_UTL_fileExists(name, true, true) == false)
    AS_UTL_mkdir(name);

  //  Open all the output files -- fail early if we can't open that many files.

  blobfiles[0]    = NULL;
  blobfileslen[0] = UINT64_MAX;
  readfiles[0]    = NULL;
  readfileslen[0] = UINT32_MAX;

  for (uint32 i=1; i<=maxPartition; i++) {
    sprintf(name,"%s/partitions/blobs.%04d", _storePath, i);

    errno = 0;
    blobfiles[i]    = fopen(name, "w");
    blobfileslen[i] = 0;

    if (errno)
      fprintf(stderr, "gkStore::gkStore_buildPartitions()-- ERROR: failed to open partition %u file '%s' for write: %s\n",
              i, name, strerror(errno)), exit(1);

    sprintf(name,"%s/partitions/reads.%04d", _storePath, i);

    errno = 0;
    readfiles[i]    = fopen(name, "w");
    readfileslen[i] = 0;

    if (errno)
      fprintf(stderr, "gkStore::gkStore_buildPartitions()-- ERROR: failed to open partition %u file '%s' for write: %s\n",
              i, name, strerror(errno)), exit(1);
  }

  //  Open the output partition map file -- we might as well fail early if we can't make it also.

  sprintf(name,"%s/partitions/map", _storePath);

  errno = 0;
  FILE *rIDmF = fopen(name, "w");
  if (errno)
    fprintf(stderr, "gkStore::gkStore_buildPartitions()-- ERROR: failed to open partition map file '%s': %s\n",
            name, strerror(errno)), exit(1);

  //  Copy the blob from the master file to the partitioned file, update pointers.

  readIDmap[0] = UINT32_MAX;    //  There isn't a zeroth read, make it bogus.

  for (uint32 fi=1; fi<=gkStore_getNumReads(); fi++) {
    uint32  pi = partitionMap[fi];

    assert(pi != 0);  //  No zeroth partition, right?

    if (pi == UINT32_MAX)
      //  Deleted reads are not assigned a partition; skip them
      continue;

    //  Make a copy of the read, then modify it for the partition, then write it to the partition.
    //  Without the copy, we'd need to update the master record too.

    gkRead  partRead = _reads[fi];  //*gkStore_getRead(fi);

    partRead.gkRead_copyDataToPartition(_blobs, blobfiles, blobfileslen, pi);

#if 1
    fprintf(stderr, "read "F_U32"="F_U32" len "F_U32" -- blob master "F_U64" -- to part "F_U32" new read id "F_U32" blob "F_U64"/"F_U64" -- at readIdx "F_U32"\n",
            fi, _reads[fi].gkRead_readID(), _reads[fi].gkRead_sequenceLength(),
            _reads[fi]._mPtr,
            pi,
            partRead.gkRead_readID(), partRead._pID, partRead._mPtr,
            readfileslen[pi]);
#endif

    AS_UTL_safeWrite(readfiles[pi], &partRead, "gkStore::gkStore_buildPartitions::read", sizeof(gkRead), 1);

    readIDmap[fi] = readfileslen[pi]++;
  }

  //  There isn't a zeroth read.

  AS_UTL_safeWrite(rIDmF, &maxPartition,  "gkStore::gkStore_buildPartitions::maxPartition", sizeof(uint32), 1);
  AS_UTL_safeWrite(rIDmF,  readfileslen,  "gkStore::gkStore_buildPartitions::readfileslen", sizeof(uint32), maxPartition + 1);
  AS_UTL_safeWrite(rIDmF,  partitionMap,  "gkStore::gkStore_buildPartitions::partitionMap", sizeof(uint32), gkStore_getNumReads() + 1);
  AS_UTL_safeWrite(rIDmF,  readIDmap,     "gkStore::gkStore_buildPartitions::readIDmap",    sizeof(uint32), gkStore_getNumReads() + 1);

  //  cleanup -- close all the files, delete storage

  fclose(rIDmF);

  for (uint32 i=1; i<=maxPartition; i++) {
    fprintf(stderr, "partition "F_U32" has "F_U32" reads\n", i, readfileslen[i]);

    errno = 0;

    fclose(blobfiles[i]);
    fclose(readfiles[i]);

    if (errno)
      fprintf(stderr, "  warning: %s\n", strerror(errno));
  }

  delete [] readIDmap;
  delete [] readfileslen;
  delete [] readfiles;
  delete [] blobfileslen;
  delete [] blobfiles;
}



void
gkStore::gkStore_delete(void) {
  char path[FILENAME_MAX];

  delete [] _libraries;
  delete [] _reads;

  _libraries = NULL;
  _reads     = NULL;

  gkStore_deletePartitions();

  sprintf(path, "%s/info",      gkStore_path());  AS_UTL_unlink(path);
  sprintf(path, "%s/libraries", gkStore_path());  AS_UTL_unlink(path);
  sprintf(path, "%s/reads",     gkStore_path());  AS_UTL_unlink(path);
  sprintf(path, "%s/blobs",     gkStore_path());  AS_UTL_unlink(path);

  AS_UTL_unlink(path);
}



void
gkStore::gkStore_deletePartitions(void) {
  char path[FILENAME_MAX];

  sprintf(path, "%s/partitions/map", gkStore_path());

  if (AS_UTL_fileExists(path, false, false) == false)
    return;

  //  How many partitions?

  FILE *F = fopen(path, "r");
  if (errno)
    fprintf(stderr, "ERROR: failed to open partition meta data '%s': %s\n", path, strerror(errno)), exit(1);

  AS_UTL_safeRead(F, &_numberOfPartitions, "gkStore_deletePartitions::numberOfPartitions", sizeof(uint32), 1);

  fclose(F);

  //  Yay!  Delete!

  AS_UTL_unlink(path);

  for (uint32 ii=0; ii<_numberOfPartitions; ii++) {
    sprintf(path, "%s/partitions/reads.%04u", gkStore_path(), ii+1);  AS_UTL_unlink(path);
    sprintf(path, "%s/partitions/blobs.%04u", gkStore_path(), ii+1);  AS_UTL_unlink(path);
  }
}






void
gkStoreStats::init(gkStore *UNUSED(gkp)) {

#if 0
  gkFragment    fr;
  gkStream     *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  numActiveFrag     = 0;
  numMatedFrag      = 0;
  readLength        = 0;
  clearLength       = 0;

  lowestID          = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  highestID         = new uint32 [gkp->gkStore_getNumLibraries() + 1];

  numActivePerLib   = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  numMatedPerLib    = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  readLengthPerLib  = new uint64 [gkp->gkStore_getNumLibraries() + 1];
  clearLengthPerLib = new uint64 [gkp->gkStore_getNumLibraries() + 1];

  for (uint32 i=0; i<gkp->gkStore_getNumLibraries() + 1; i++) {
    lowestID[i]          = 0;
    highestID[i]         = 0;

    numActivePerLib[i]   = 0;
    numMatedPerLib[i]    = 0;
    readLengthPerLib[i]  = 0;
    clearLengthPerLib[i] = 0;
  }

  while (fs->next(&fr)) {
    uint32     lid = fr.gkFragment_getLibraryID();
    uint32     rid = fr.gkFragment_getReadID();

    if (lowestID[lid] == 0) {
      lowestID[lid]  = rid;
      highestID[lid] = rid;
    }
    if (highestID[lid] < rid) {
      highestID[lid] = rid;
    }

    numActiveFrag++;
    numActivePerLib[lid]++;

    readLength             += fr.gkFragment_getSequenceLength();
    readLengthPerLib[lid]  += fr.gkFragment_getSequenceLength();

    clearLength            += fr.gkFragment_getClearRegionLength();
    clearLengthPerLib[lid] += fr.gkFragment_getClearRegionLength();
  }

  delete fs;
#endif

}

