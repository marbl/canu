#include "util++.H"

#include "tapperTag.H"
#include "tapperResult.H"
#include "tapperAlignment.H"
#include "tapperHit.H"
#include "tapperGlobalData.H"
#include "tapperThreadData.H"
#include "tapperComputation.H"



//  Reads a tapperAlignmentFile, converts all the alignments to
//  tapperAlignments (loses mate pair information), and sorts by
//  position on the reference.

//  There are (at least) two ways to sort.  Merge sort or bucket sort.
//
//  Bucket sort is a little easier, but, without knowing the length of
//  the reference sequences, we cannot map seq,pos to a bucket.  We
//  also have no memory guarantee; it is possible to have a bucket get
//  too big.
//
//  Merge sort is more difficult, because of the merge.  We have a
//  memory size guarantee though.


uint32
saveFrag(tapperAlignment *ali, uint32 aliLen, tapperResult *res, uint32 fragLen, tapperResultFragment *frag) {

  for (uint32 i=0; i<fragLen; i++) {
    tapperResultFragment *f = frag + i;

    //  At least one is true, and at least one is false ==> exactly
    //  one is true.

    if ((f->_qual._tag1valid == 0) && (f->_qual._tag2valid == 0))
      fprintf(stderr, "error\n");

    assert((f->_qual._tag1valid == 1) || (f->_qual._tag2valid == 1));
    assert((f->_qual._tag1valid == 0) || (f->_qual._tag2valid == 0));

    if (f->_qual._tag1valid) {
      memset(ali + aliLen, 0, sizeof(tapperAlignment));

      ali[aliLen]._tagid              = res->idx._tag1id;
      ali[aliLen]._seq                = f->_seq;
      ali[aliLen]._pos                = f->_pos;
      ali[aliLen]._basesMismatch      = f->_qual._tag1basesMismatch;
      ali[aliLen]._colorMismatch      = f->_qual._tag1colorMismatch;
      ali[aliLen]._colorInconsistent  = f->_qual._tag1colorInconsistent;
      ali[aliLen]._rev                = f->_qual._tag1rev;

      ali[aliLen]._diffSize           = f->_qual._diffSize;

      memcpy(ali[aliLen]._colorDiffs,
             f->_qual._tag1colorDiffs,
             sizeof(uint8) * MAX_COLOR_MISMATCH_MAPPED);

      aliLen++;
    }

    if (f->_qual._tag2valid) {
      memset(ali + aliLen, 0, sizeof(tapperAlignment));

      ali[aliLen]._tagid              = res->idx._tag2id;
      ali[aliLen]._seq                = f->_seq;
      ali[aliLen]._pos                = f->_pos;
      ali[aliLen]._basesMismatch      = f->_qual._tag2basesMismatch;
      ali[aliLen]._colorMismatch      = f->_qual._tag2colorMismatch;
      ali[aliLen]._colorInconsistent  = f->_qual._tag2colorInconsistent;
      ali[aliLen]._rev                = f->_qual._tag2rev;

      ali[aliLen]._diffSize           = f->_qual._diffSize;

      memcpy(ali[aliLen]._colorDiffs,
             f->_qual._tag2colorDiffs,
             sizeof(uint8) * MAX_COLOR_MISMATCH_MAPPED);

      aliLen++;
    }
  }

  return(aliLen);
}



uint32
saveMate(tapperAlignment *ali, uint32 aliLen, tapperResult *res) {

  for (uint32 i=0; i<res->idx._numMated; i++) {
    tapperResultMated *m = res->mate + i;

    memset(ali + aliLen, 0, sizeof(tapperAlignment));

    ali[aliLen]._tagid              = res->idx._tag1id;
    ali[aliLen]._seq                = m->_seq;
    ali[aliLen]._pos                = m->_pos1;
    ali[aliLen]._basesMismatch      = m->_qual._tag1basesMismatch;
    ali[aliLen]._colorMismatch      = m->_qual._tag1colorMismatch;
    ali[aliLen]._colorInconsistent  = m->_qual._tag1colorInconsistent;
    ali[aliLen]._rev                = m->_qual._tag1rev;

    ali[aliLen]._diffSize           = m->_qual._diffSize;

    memcpy(ali[aliLen]._colorDiffs,
           m->_qual._tag1colorDiffs,
           sizeof(uint8) * MAX_COLOR_MISMATCH_MAPPED);

    aliLen++;

    memset(ali + aliLen, 0, sizeof(tapperAlignment));

    ali[aliLen]._tagid              = res->idx._tag2id;
    ali[aliLen]._seq                = m->_seq;
    ali[aliLen]._pos                = m->_pos2;
    ali[aliLen]._basesMismatch      = m->_qual._tag2basesMismatch;
    ali[aliLen]._colorMismatch      = m->_qual._tag2colorMismatch;
    ali[aliLen]._colorInconsistent  = m->_qual._tag2colorInconsistent;
    ali[aliLen]._rev                = m->_qual._tag2rev;

    ali[aliLen]._diffSize           = m->_qual._diffSize;

    memcpy(ali[aliLen]._colorDiffs,
           m->_qual._tag2colorDiffs,
           sizeof(uint8) * MAX_COLOR_MISMATCH_MAPPED);

    aliLen++;
  }

  return(aliLen);
}




uint32
sortAndDump(tapperAlignment *ali, uint32 aliLen, char *outputName, uint32 &outputIndex) {
  char     filename[FILENAME_MAX];

  if (aliLen == 0)
    return(0);

  tapperAlignmentPositionCompare pc;
  std::sort(ali, ali + aliLen, pc);

  sprintf(filename, "%s."uint32FMTW(03)".tapperAlignment", outputName, outputIndex);

  fprintf(stderr, "Writing "uint32FMT" sorted alignments to '%s'\n", aliLen, filename);

  recordFile  *out = new recordFile(filename, 0, sizeof(tapperAlignment), 'w');
  out->putRecord(ali, aliLen);
  delete out;

  outputIndex++;

  return(0);
}



int
main(int argc, char **argv) {
  char     *outputName  = 0L;
  uint32    outputIndex = 0;
  uint32    inputsLen   = 0;
  char     *inputs[8192];
  uint64    memoryLimit = 1024 * 1024 * 1024;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-memory", 2) == 0) {
      memoryLimit = strtouint64(argv[++arg], 0L) * 1024 * 1024;

    } else if (strncmp(argv[arg], "-output", 2) == 0) {
      outputName = argv[++arg];

    } else {
      if (tapperResultFile::validResultFile(argv[arg]) == false) {
        fprintf(stderr, "Didn't find tapperResultFile '%s'\n", argv[arg]);
        err++;
      } else {
        inputs[inputsLen++] = argv[arg];
      }
    }

    arg++;
  }
  if ((err) || (inputsLen == 0) || (outputName == 0L)) {
    fprintf(stderr, "usage: %s [-memory X (MB)] -output prefix input ....\n", argv[0]);
    exit(1);
  }


  {
    uint32              aliMax = memoryLimit / sizeof(tapperAlignment);
    uint32              aliLen = 0;
    tapperAlignment    *ali    = new tapperAlignment [aliMax];

    fprintf(stderr, "Can fit "uint32FMT" alignments into "uint64FMT" bytes memory; "uint32FMT" bytes each.\n",
            aliMax, memoryLimit, (uint32)sizeof(tapperAlignment));

    speedCounter        S(" %10.0f results (%8.0f results/sec)\r", 1, 100000, true);

    for (uint32 inputsIdx=0; inputsIdx<inputsLen; inputsIdx++) {
      tapperResultFile   *inp = new tapperResultFile(inputs[inputsIdx], 'r');
      tapperResult       *res = new tapperResult;

      while (inp->read(res)) {

        //  Sort and dump if the next result has too many alignments.
        //
        if (aliMax < aliLen + (res->idx._numFrag +
                               res->idx._numFragSingleton +
                               res->idx._numFragTangled +
                               res->idx._numMated * 2)) {
          aliLen = sortAndDump(ali, aliLen, outputName, outputIndex);
        }

        aliLen = saveFrag(ali, aliLen, res, res->idx._numFrag,          res->frag);
        aliLen = saveFrag(ali, aliLen, res, res->idx._numFragSingleton, res->sing);
        aliLen = saveFrag(ali, aliLen, res, res->idx._numFragTangled,   res->tali);
        aliLen = saveMate(ali, aliLen, res);

        S.tick();
      }
      S.finish();

      delete inp;
      delete res;
    }

    aliLen = sortAndDump(ali, aliLen, outputName, outputIndex);

    delete [] ali;
  }

  //
  //  Now the merge.
  //

  {
    char                filename[FILENAME_MAX];

    tapperAlignment    *ali = new tapperAlignment   [outputIndex];
    recordFile        **inp = new recordFile      * [outputIndex];
    recordFile         *out = 0L;

    bool                stillMore = true;
    uint32              minidx = 0;

    tapperAlignmentPositionCompare lessthan;

    for (uint32 x=0; x<outputIndex; x++) {
      sprintf(filename, "%s."uint32FMTW(03)".tapperAlignment", outputName, x);
      inp[x] = new recordFile(filename, 0, sizeof(tapperAlignment), 'r');

      inp[x]->getRecord(ali + x);
    }

    sprintf(filename, "%s.tapperAlignment", outputName);
    out = new recordFile(filename, 0, sizeof(tapperAlignment), 'w');


    while (stillMore) {

      //  Compare all against the current default minidx, pick the
      //  smallest alignment currently loaded.
      for (uint32 x=0; x<outputIndex; x++)
        if ((x != minidx) && (inp[x] != 0L) && (lessthan(ali[x], ali[minidx])))
          minidx = x;

      //  Dump it.
      out->putRecord(ali + minidx);

      //  Read the next record.  If no next record, close the file,
      //  and pick a new default minidx
      if (inp[minidx]->getRecord(ali + minidx) == 0) {
        delete inp[minidx];
        inp[minidx] = 0L;

        stillMore = false;

        for (uint32 x=0; x<outputIndex; x++)
          if (inp[x] != 0L) {
            minidx    = x;
            stillMore = true;
          }
      }
    }

    delete    out;

    for (uint32 x=0; x<outputIndex; x++) {
      assert(inp[x] == 0L);

      sprintf(filename, "%s."uint32FMTW(03)".tapperAlignment", outputName, x);
      unlink(filename);
    }

    delete [] inp;
    delete [] ali;
  }

  exit(0);
}
