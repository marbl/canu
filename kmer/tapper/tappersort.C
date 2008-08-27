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

int
main(int argc, char **argv) {
  char     *outputName = 0L;
  char     *resultName = 0L;

  u64bit    memoryLimit = 1024 * 1024 * 1024;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-memory", 2) == 0) {
      memoryLimit = strtou64bit(argv[++arg], 0L) * 1024 * 1024;

    } else if (strncmp(argv[arg], "-output", 2) == 0) {
      outputName = argv[++arg];

    } else if (resultName == 0L) {
      resultName = argv[arg];

    } else {
      err++;
    }

    arg++;
  }
  if ((err) || (resultName == 0L) || (outputName == 0L)) {
    fprintf(stderr, "usage: %s [-memory X (MB)] -output prefix input\n", argv[0]);
    exit(1);
  }

  tapperResultFile   *inp = new tapperResultFile(resultName, 'r');
  tapperResult       *res = new tapperResult;

  u32bit              aliMax = memoryLimit / sizeof(tapperAlignment);
  u32bit              aliLen = 0;
  tapperAlignment    *ali    = new tapperAlignment [aliMax];

  speedCounter        S(" %8f alignments (%8.5f alignments/sec)\r", 1, 100000, true);

  while (inp->read(res)) {

    //  Sort and dump if the next result has too many alignments.
    //
    if (aliMax < aliLen + (res->idx._numFragment +
                           res->idx._numSingleton +
                           res->idx._numMated * 2)) {
    }



    for (u32bit i=0; i<res->idx._numFragment; i++) {
      tapperResultFragment *f = res->frag + i;

      assert(f->_qual._tag1valid == true);
      assert(f->_qual._tag2valid == false);

      ali[aliLen]._tagid              = res->idx._tag1id;
      ali[aliLen]._seq                = f->_seq;
      ali[aliLen]._pos                = f->_pos;
      ali[aliLen]._pad                = 0;
      ali[aliLen]._basesMismatch      = f->_qual._tag1basesMismatch;
      ali[aliLen]._colorMismatch      = f->_qual._tag1colorMismatch;
      ali[aliLen]._colorInconsistent  = f->_qual._tag1colorInconsistent;
      ali[aliLen]._rev                = f->_qual._tag1rev;

      ali[aliLen]._diffSize           = f->_qual._diffSize;

      memcpy(ali[aliLen]._colorDiffs,
             f->_qual._tag1colorDiffs,
             sizeof(u8bit) * MAX_COLOR_MISMATCH_MAPPED);

      aliLen++;
      S.tick();
    }

    for (u32bit i=0; i<res->idx._numSingleton; i++) {
      tapperResultSingleton *s = res->sing + i;

      //  At least one is true, and at least one is false ==> exactly
      //  one is true.

      assert((s->_qual._tag1valid == true)  || (s->_qual._tag2valid == true));
      assert((s->_qual._tag1valid == false) || (s->_qual._tag2valid == false));

      if (s->_qual._tag2valid) {
        ali[aliLen]._tagid              = res->idx._tag1id;
        ali[aliLen]._seq                = s->_seq;
        ali[aliLen]._pos                = s->_pos;
        ali[aliLen]._pad                = 0;
        ali[aliLen]._basesMismatch      = s->_qual._tag1basesMismatch;
        ali[aliLen]._colorMismatch      = s->_qual._tag1colorMismatch;
        ali[aliLen]._colorInconsistent  = s->_qual._tag1colorInconsistent;
        ali[aliLen]._rev                = s->_qual._tag1rev;

        ali[aliLen]._diffSize           = s->_qual._diffSize;

        memcpy(ali[aliLen]._colorDiffs,
               s->_qual._tag1colorDiffs,
               sizeof(u8bit) * MAX_COLOR_MISMATCH_MAPPED);

        aliLen++;
        S.tick();
      }

      if (s->_qual._tag2valid) {
        ali[aliLen]._tagid              = res->idx._tag2id;
        ali[aliLen]._seq                = s->_seq;
        ali[aliLen]._pos                = s->_pos;
        ali[aliLen]._pad                = 0;
        ali[aliLen]._basesMismatch      = s->_qual._tag2basesMismatch;
        ali[aliLen]._colorMismatch      = s->_qual._tag2colorMismatch;
        ali[aliLen]._colorInconsistent  = s->_qual._tag2colorInconsistent;
        ali[aliLen]._rev                = s->_qual._tag2rev;

        ali[aliLen]._diffSize           = s->_qual._diffSize;

        memcpy(ali[aliLen]._colorDiffs,
               s->_qual._tag2colorDiffs,
               sizeof(u8bit) * MAX_COLOR_MISMATCH_MAPPED);

        aliLen++;
        S.tick();
      }
    }

    for (u32bit i=0; i<res->idx._numMated; i++) {
      tapperResultMated *m = res->mate + i;

      ali[aliLen]._tagid              = res->idx._tag1id;
      ali[aliLen]._seq                = m->_seq;
      ali[aliLen]._pos                = m->_pos1;
      ali[aliLen]._pad                = 0;
      ali[aliLen]._basesMismatch      = m->_qual._tag1basesMismatch;
      ali[aliLen]._colorMismatch      = m->_qual._tag1colorMismatch;
      ali[aliLen]._colorInconsistent  = m->_qual._tag1colorInconsistent;
      ali[aliLen]._rev                = m->_qual._tag1rev;

      ali[aliLen]._diffSize           = m->_qual._diffSize;

      memcpy(ali[aliLen]._colorDiffs,
             m->_qual._tag1colorDiffs,
             sizeof(u8bit) * MAX_COLOR_MISMATCH_MAPPED);

      aliLen++;
      S.tick();


      ali[aliLen]._tagid              = res->idx._tag2id;
      ali[aliLen]._seq                = m->_seq;
      ali[aliLen]._pos                = m->_pos2;
      ali[aliLen]._pad                = 0;
      ali[aliLen]._basesMismatch      = m->_qual._tag2basesMismatch;
      ali[aliLen]._colorMismatch      = m->_qual._tag2colorMismatch;
      ali[aliLen]._colorInconsistent  = m->_qual._tag2colorInconsistent;
      ali[aliLen]._rev                = m->_qual._tag2rev;

      ali[aliLen]._diffSize           = m->_qual._diffSize;

      memcpy(ali[aliLen]._colorDiffs,
             m->_qual._tag2colorDiffs,
             sizeof(u8bit) * MAX_COLOR_MISMATCH_MAPPED);

      aliLen++;
      S.tick();
    }

    //for (u32bit i=0; i<res->idx._numTangled; i++) {
    //  res->tang[i].print(stdout, &res->idx);
    //}
  }

  S.finish();

  delete inp;
  delete res;

  exit(0);
}
