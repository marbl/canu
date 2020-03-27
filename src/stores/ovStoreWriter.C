
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "ovStore.H"


////////////////////////////////////////
//
//  SEQUENTIAL STORE - only two functions.
//

ovStoreWriter::ovStoreWriter(const char *path, sqStore *seq) {
  char name[FILENAME_MAX+1];

  memset(_storePath, 0, FILENAME_MAX);
  strncpy(_storePath, path, FILENAME_MAX);

  //  Fail if this is a valid ovStore.
#warning not failing if a valid store
  //if (_info.test(_storePath) == true)
  //  fprintf(stderr, "ERROR:  '%s' is a valid ovStore; cannot create a new one.\n", _storePath), exit(1);

  //  Create the new store

  AS_UTL_mkdir(_storePath);

  _info.clear(seq->sqStore_lastReadID());
  //_info.save(_storePath);   Used to save this as a sentinel, but now fails asserts I like

  _seq       = seq;

  _index     = new ovStoreOfft [_info.maxID() + 1];

  _bof       = NULL;   //  Open the file on the first overlap.
  _bofSlice  = 1;      //  Constant, never changes.
  _bofPiece  = 1;      //  Incremented whenever a file is closed.

  _histogram = new ovStoreHistogram(_seq);  //  Only used for merging in results from output files.
}



ovStoreWriter::~ovStoreWriter() {

  //  Write the index

  AS_UTL_saveFile(_storePath, '/', "index", _index, _info.maxID()+1);

  delete [] _index;

  //  Update our copy of the histogram from the last open file, and close it.

  if (_bof) {
    _histogram->mergeHistogram(_bof->getHistogram());
    _bof->removeHistogram();
  }

  delete _bof;

  //  Save the histogram data.

  _histogram->saveHistogram(_storePath);

  delete _histogram;

  //  Update the on-disk info with the results and real magic number

  _info.save(_storePath);

  //  Report a nice success message.

  fprintf(stderr, "Created ovStore '%s' with " F_U64 " overlaps for reads from " F_U32 " to " F_U32 ".\n",
          _storePath, _info.numOverlaps(), _info.bgnID(), _info.endID());
}



void
ovStoreWriter::writeOverlap(ovOverlap *overlap) {

  //  Close the current output file if it's too big.
  //    The current output file must exist.
  //    The current output file must be too big.
  //    The current output file must NOT contain overlaps for the same read this overlap does.
  //

  //fprintf(stderr, "writeOverlap() maxID %u endID %u thisID %u\n", _info.maxID(), _info.endID(), overlap->a_iid);

  //  If an output file and it's too big (and we're not in the middle of a set of overlaps
  //  for some read), make a new output file.  Save the histogram from the previous
  //  file in our master list, then delete it so it isn't saved to disk.

  if ((_bof != NULL) &&
      (_bof->fileTooBig() == true) &&
      (_info.endID() < overlap->a_iid)) {
    _histogram->mergeHistogram(_bof->getHistogram());

    _bof->removeHistogram();

    delete _bof;

    _bof = NULL;
    _bofPiece++;
  }

  //  Open a new output file if there isn't one.

  if (_bof == NULL)
    _bof = new ovFile(_seq, _storePath, _bofSlice, _bofPiece, ovFileNormalWrite);

  //  Make sure the overlaps are sorted, and add the overlap to the info file.

  if (overlap->a_iid > _info.maxID()) {
    assert(0);
  }

  //  Add the overlap to the index and info.

  _index[overlap->a_iid].addOverlap(_bofSlice, _bofPiece, _bof->filePosition(), _info.numOverlaps());

  _info.addOverlaps(overlap->a_iid, 1);

  //  Write the overlap.

  _bof->writeOverlap(overlap);
}




////////////////////////////////////////
//
//  PARALLEL STORE - many functions, all the rest.
//

ovStoreSliceWriter::ovStoreSliceWriter(const char *path,
                                       sqStore    *seq,
                                       uint32      sliceNum,
                                       uint32      numSlices,
                                       uint32      numBuckets) {

  memset(_storePath, 0, FILENAME_MAX);
  strncpy(_storePath, path, FILENAME_MAX);

  _seq                 = seq;

  _sliceNum            = sliceNum;
  _pieceNum            = 1;
  _numSlices           = numSlices;
  _numBuckets          = numBuckets;
};



ovStoreSliceWriter::~ovStoreSliceWriter() {

  //  Report a nice success message.

  //fprintf(stderr, "Created ovStore '%s' with " F_U64 " overlaps for reads from " F_U32 " to " F_U32 ".\n",
  //        _storePath, _info.numOverlaps(), _info.bgnID(), _info.endID());
}




uint64
ovStoreSliceWriter::loadBucketSizes(uint64 *bucketSizes) {
  char      name[FILENAME_MAX+1];

  uint64   *sliceSizes = new uint64 [_numSlices + 1];  //  For each overlap job, number of overlaps per bucket
  uint64    totOvl     = 0;

  for (uint32 i=0; i<=_numBuckets; i++) {
    bucketSizes[i] = 0;

    snprintf(name, FILENAME_MAX, "%s/bucket%04u/slice%04u", _storePath, i, _sliceNum);

    //  If no file, there are no overlaps, so nothing to load.

    if (fileExists(name) == false)
      continue;

    //  Load the slice sizes, and save the number of overlaps in this slice.

    snprintf(name, FILENAME_MAX, "%s/bucket%04u/sliceSizes", _storePath, i);
    AS_UTL_loadFile(name, sliceSizes, _numSlices + 1);  //  Checks that all data is loaded, too.

    fprintf(stderr, "  found %10" F_U64P " overlaps in '%s'.\n", sliceSizes[_sliceNum], name);

    bucketSizes[i] = sliceSizes[_sliceNum];
    totOvl        += sliceSizes[_sliceNum];
  }

  delete [] sliceSizes;

  return(totOvl);
}



void
ovStoreSliceWriter::loadOverlapsFromBucket(uint32 bucket, uint64 expectedLen, ovOverlap *ovls, uint64& ovlsLen) {
  char name[FILENAME_MAX+1];

  if (expectedLen == 0)
    return;

  snprintf(name, FILENAME_MAX, "%s/bucket%04u/slice%04u", _storePath, bucket, _sliceNum);

  if (fileExists(name) == false)
    fprintf(stderr, "ERROR: " F_U64 " overlaps claim to exist in bucket '%s', but file not found.\n",
            expectedLen, name), exit(1);

  fprintf(stderr, "  loading %10" F_U64P " overlaps from '%s'.\n", expectedLen, name);

  ovFile   *bof    = new ovFile(_seq, name, ovFileFull);
  uint64    before = ovlsLen;

  while (bof->readOverlap(ovls + ovlsLen))
    ovlsLen++;

  delete bof;

  if (ovlsLen - before != expectedLen)
    fprintf(stderr, "ERROR: expected " F_U64 " overlaps, found " F_U64 " overlaps.\n",
            expectedLen, ovlsLen - before), exit(1);
}



void
ovStoreSliceWriter::writeOverlaps(ovOverlap  *ovls,
                                  uint64      ovlsLen) {
  ovStoreInfo    info(_seq->sqStore_lastReadID());

  //  Probably wouldn't be too hard to make this take all overlaps for one read.
  //  But would need to track the open files in the class, not only in this function.
  assert(info.numOverlaps() == 0);

  //  Check that overlaps are sorted.

  uint64  nUnsorted = 0;

  for (uint64 oo=1; oo<ovlsLen; oo++)
    if (ovls[oo-1].a_iid > ovls[oo].a_iid)
      nUnsorted++;

  if (nUnsorted > 0) {
    fprintf(stderr, "ERROR: Overlaps aren't sorted.\n");
    exit(1);
  }

  //  Create the index and overlaps files

  ovStoreOfft  *index     = new ovStoreOfft [_seq->sqStore_lastReadID() + 1];
  ovFile       *olapFile  = new ovFile(_seq, _storePath, _sliceNum, _pieceNum, ovFileNormalWrite);

  //  Dump the overlaps

  for (uint64 oo=0; oo<ovlsLen; oo++ ) {

    //  If this overlap is for a new read, and we've written too many overlaps
    //  to the current piece, start a new piece.

    if ((olapFile->fileTooBig() == true) &&
        (ovls[oo].a_iid          > info.endID())) {
      delete olapFile;

      _pieceNum++;

      olapFile  = new ovFile(_seq, _storePath, _sliceNum, _pieceNum, ovFileNormalWrite);
    }

    //  Add the overlap to the index.

    index[ovls[oo].a_iid].addOverlap(_sliceNum, _pieceNum, olapFile->filePosition(), oo);

    //  Add the overlap to the file.

    olapFile->writeOverlap(ovls + oo);

    //  Add the overlap to the info

    info.addOverlaps(ovls[oo].a_iid, 1);
  }

  //  Close the output file, write the index, write the info.

  delete    olapFile;

  char indexName[FILENAME_MAX+1];
  snprintf(indexName, FILENAME_MAX, "%s/%04u.index", _storePath, _sliceNum);
  AS_UTL_saveFile(indexName, index, info.maxID()+1);

  delete [] index;

  info.save(_storePath, _sliceNum, true);

  //  And done.

  fprintf(stderr, "  created '%s/%04u' with " F_U64 " overlaps for reads " F_U32 " to " F_U32 ".\n",
          _storePath, _sliceNum, info.numOverlaps(), info.bgnID(), info.endID());
}





void
ovStoreSliceWriter::mergeInfoFiles(void) {
  char           indexName[FILENAME_MAX+1];

  ovStoreInfo    infopiece[_numSlices + 1];

  //  Load the info files.  We need to load at least one to figure out maxID.

  fprintf(stderr, " - Loading slice information.\n");
  fprintf(stderr, " -\n");
  fprintf(stderr, " - slice     bgnID     endID    maxID   numOlaps\n");
  fprintf(stderr, " - ----- --------- --------- -------- ----------\n");

  for (uint32 ss=1; ss<=_numSlices; ss++) {
    infopiece[ss].load(_storePath, ss, true);

    fprintf(stderr, " - %5" F_U32P " %9" F_U32P " %9" F_U32P " %8" F_U32P " %10" F_U64P "\n",
            ss,
            infopiece[ss].bgnID(),
            infopiece[ss].endID(),
            infopiece[ss].maxID(),
            infopiece[ss].numOverlaps());
  }

  fprintf(stderr, " - ----- --------- --------- -------- ----------\n");

  //  Set us up and allocate some space.

  ovStoreInfo    info(infopiece[1].maxID());

  ovStoreOfft   *indexpiece = new ovStoreOfft [infopiece[1].maxID() + 1];
  ovStoreOfft   *index      = new ovStoreOfft [infopiece[1].maxID() + 1];

  //  Merge in indexes.

  fprintf(stderr, " -\n");
  fprintf(stderr, " - Merging indexes.\n");
  fprintf(stderr, " -\n");

  for (uint32 ss=1; ss<=_numSlices; ss++) {

    //  Load the index for this piece.

    snprintf(indexName, FILENAME_MAX, "%s/%04u.index", _storePath, ss);
    AS_UTL_loadFile(indexName, indexpiece, info.maxID()+1);

    //  Copy index elements to the master index, updating overlapID.

    for (uint32 ii=infopiece[ss].bgnID(); ii<=infopiece[ss].endID(); ii++) {
      index[ii]             = indexpiece[ii];
      index[ii]._overlapID += info.numOverlaps();
    }

    //  Update the master info.

    info.addOverlaps(infopiece[ss].bgnID(), 0);
    info.addOverlaps(infopiece[ss].endID(), infopiece[ss].numOverlaps());

    //  Keep users entertained.

    //fprintf(stderr, " -  Now finished with fragments " F_U32 " to " F_U32 " -- " F_U64 " overlaps.\n",
    //        info.bgnID(), info.endID(), info.numOverlaps());
  }

  //  Dump the new store info file

  info.save(_storePath, _numSlices);

  //  And the new index

  AS_UTL_saveFile(_storePath, '/', "index", index, info.maxID()+1);

  //  Cleanup and done!

  delete [] indexpiece;
  delete [] index;

  fprintf(stderr, " - Finished.  " F_U32 " reads with " F_U64 " overlaps.\n",
          info.endID(), info.numOverlaps());
  fprintf(stderr, " -\n");
}



void
ovStoreSliceWriter::mergeHistogram(void) {
  char               dataname[FILENAME_MAX+1] = {0};

  fprintf(stderr, " - Merging histograms.\n");
  fprintf(stderr, " -\n");
  fprintf(stderr, " -\n");
  fprintf(stderr, " - slice piece     bgnID     endID\n");
  fprintf(stderr, " - ----- ----- --------- ---------\n");

  ovStoreHistogram  *merged = new ovStoreHistogram(_seq);

  for (uint32 ss=1; ss <= _numSlices; ss++) {
    for (uint32 pp=1; pp < 1000; pp++) {
      ovStoreHistogram  piece(ovFile::createDataName(dataname, _storePath, ss, pp));

      if (piece.overlapScoresLastID() > 0) {
        fprintf(stderr, " - %5u %5u %9u %9u\n",
                ss, pp, piece.overlapScoresBaseID(), piece.overlapScoresLastID());

        merged->mergeHistogram(&piece);
      } else {
        break;
      }
    }
  }

  merged->saveHistogram(_storePath);

  fprintf(stderr, " - ----- ----- --------- ---------\n");

  delete merged;
}



void
ovStoreSliceWriter::removeOverlapSlice(void) {
  char name[FILENAME_MAX+1];

  for (uint32 bb=0; bb<=_numBuckets; bb++) {
    snprintf(name, FILENAME_MAX, "%s/bucket%04u/slice%04u", _storePath, bb, _sliceNum);
    AS_UTL_unlink(name);
  }
}



void
ovStoreSliceWriter::checkSortingIsComplete(void) {
  char    nameF[FILENAME_MAX+1];
  char    nameI[FILENAME_MAX+1];

  fprintf(stderr, " - Checking that sorting is complete (every slice should have an info file).\n");
  fprintf(stderr, " -\n");
  fprintf(stderr, "\n");

  uint32  failedJobs = 0;

  for (uint32 i=1; i<=_numSlices; i++) {
    snprintf(nameF, FILENAME_MAX, "%s/%04u.info",  _storePath, i);
    snprintf(nameI, FILENAME_MAX, "%s/%04u.index", _storePath, i);

    bool existF = fileExists(nameF);
    bool existI = fileExists(nameI);

    if (existF && existI)
      continue;

    failedJobs++;

    if (existF == false)    fprintf(stderr, "ERROR: Segment " F_U32 " info  not present (%s)\n", i, nameF);
    if (existI == false)    fprintf(stderr, "ERROR: Segment " F_U32 " index not present (%s)\n", i, nameI);
  }

  if (failedJobs > 0)
    fprintf(stderr, "ERROR: " F_U32 " segments, out of " F_U32 ", failed.\n", _numSlices, failedJobs), exit(1);
}



void
ovStoreSliceWriter::removeAllIntermediateFiles(void) {
  char name[FILENAME_MAX+1];
  char nomo[FILENAME_MAX+1];   //  Esperanto, in case you were wondering.

  //  Remove sorting (slices) intermediates.

  fprintf(stderr, " - Removing intermediate files.\n");
  fprintf(stderr, " -\n");

  for (uint32 ss=1; ss <= _numSlices; ss++) {
    snprintf(name, FILENAME_MAX, "%s/%04u.index", _storePath, ss);
    AS_UTL_unlink(name);

    snprintf(name, FILENAME_MAX, "%s/%04u.info",  _storePath, ss);
    AS_UTL_unlink(name);

    for (uint32 pp=1; pp < 1000; pp++) {
      ovFile::createDataName(name, _storePath, ss, pp);
      ovStoreHistogram::createDataName(nomo, name);

      if (fileExists(nomo) == false)
        break;

      AS_UTL_unlink(nomo);
    }
  }

  //  Remove buckets.

  for (uint32 bb=1; bb <= _numBuckets; bb++) {
    snprintf(name, FILENAME_MAX, "%s/bucket%04u/sliceSizes", _storePath, bb);
    AS_UTL_unlink(name);

    for (uint32 ss=1; ss <= _numSlices; ss++) {
      snprintf(name, FILENAME_MAX, "%s/bucket%04u/slice%04u", _storePath, bb, ss);
      AS_UTL_unlink(name);
    }

    snprintf(name, FILENAME_MAX, "%s/bucket%04u",            _storePath, bb);
    AS_UTL_rmdir(name);
  }
}
