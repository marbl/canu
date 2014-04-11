#include "util++.H"

#include "tapperTag.H"
#include "tapperResult.H"
#include "tapperAlignment.H"
#include "tapperHit.H"
#include "tapperGlobalData.H"
#include "tapperThreadData.H"
#include "tapperComputation.H"



class alignmentList {
public:
  alignmentList(recordFile *inp) {
    alignsMax      = 16;
    aligns         = new tapperAlignment * [alignsMax];
    alignsLen      = new uint32            [alignsMax];
    alignsPerBlock = 16384;
    alignsInp      = inp;

    for (uint32 i=0; i<alignsMax; i++) {
      aligns[i] = new tapperAlignment [alignsPerBlock];
      alignsLen[i] = alignsInp->getRecord(aligns[i], alignsPerBlock);
      fprintf(stderr, "block "uint32FMT" has "uint32FMT" things.\n", i, alignsLen[i]);
    }
  };

  ~alignmentList() {
    for (uint32 i=0; i<alignsMax; i++)
      delete [] aligns[i];

    delete [] aligns;
  };

  //  If the last element in the first block is below the specified
  //  seq,pos, we can dump all those alignments and get more.
  //
  void  trimBeforeSeqPos(uint32 seq, uint32 pos) {

  trimBeforeSeqPosAgain:
    if (alignsLen[0] == 0)
      return;

    if ((aligns[0][alignsLen[0]-1]._seq <= seq) &&
        (aligns[0][alignsLen[0]-1]._pos <  pos)) {
      tapperAlignment *save = aligns[0];

      fprintf(stderr, "block[0] - seq "uint32FMT" pos "uint32FMT"\n",
              aligns[0][alignsLen[0]-1]._seq,
              aligns[0][alignsLen[0]-1]._pos);

      for (uint32 i=1; i<alignsMax; i++) {
        aligns[i-1]    = aligns[i];
        alignsLen[i-1] = alignsLen[i];
      }

      aligns[alignsMax-1] = save;

      alignsLen[alignsMax-1] = alignsInp->getRecord(aligns[alignsMax-1], alignsPerBlock);

      fprintf(stderr, "block "uint32FMT" has "uint32FMT" things.\n", alignsMax-1, alignsLen[alignsMax-1]);

      goto trimBeforeSeqPosAgain;
    }
  };

  tapperAlignment *operator[](uint32 x) {
    uint32 block = x / alignsPerBlock;
    uint32 piece = x % alignsPerBlock;

    if (piece < alignsLen[block])
      return(aligns[block] + piece);

    return(0L);
  };

  bool empty(void) {
    return(alignsLen[0] == 0);
  };

private:
  uint32            alignsMax;

  tapperAlignment **aligns;
  uint32           *alignsLen;

  uint32            alignsPerBlock;

  recordFile       *alignsInp;
};





int
main(int argc, char **argv) {
  char     *outputName  = 0L;
  char     *inputName   = 0L;

  uint64    memoryLimit = 1024 * 1024 * 1024;

  {
  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-memory", 2) == 0) {
      memoryLimit = strtouint64(argv[++arg], 0L) * 1024 * 1024;

    } else if (strncmp(argv[arg], "-output", 2) == 0) {
      outputName = argv[++arg];

    } else if (strncmp(argv[arg], "-input", 2) == 0) {
      inputName = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if ((err) || (inputName == 0) || (outputName == 0L)) {
    fprintf(stderr, "usage: %s [-memory X (MB)] -output prefix -input inp.tapperAlignment\n", argv[0]);
    exit(1);
  }
  }

  recordFile       *inp = new recordFile(inputName, 0, sizeof(tapperAlignment), 'r');
  alignmentList     all(inp);

  uint32            winSz = 200;
  uint32            winLo = 0;
  uint32            winHi = winLo + winSz;

  uint32            linesMax = 1024;

  char              lines[1024][256];
  uint32            lineLen[1024];

  uint16            id[4];

  while (all.empty() == false) {
    memset(lines, ' ', sizeof(char) * linesMax * 256);

    for (uint32 i=0; i<linesMax; i++)
      lineLen[i] = 0;

    for (uint32 a=0; (all[a] != 0L) && (all[a]->_pos < winHi); a++) {
      tapperAlignment *rec = all[a];

      //  XXX  we lose reads that wrap into our region

      if (winLo < rec->_pos) {
        for (uint32 l=0; l<linesMax; l++) {
        if (lineLen[l] < rec->_pos - winLo) {

          //fprintf(stdout, "at l="uint32FMT" x="uint32FMT" len="uint32FMT"\n", l, rec->_pos - winLo, lineLen[l]);

#warning need the real read size here

          for (uint32 x=rec->_pos - winLo; x<rec->_pos - winLo + 25; x++)
            lines[l][x] = '.';

          //  Needed so we can disable ID printing.
          lines[l][rec->_pos - winLo + 25] = 0;

#undef WITH_IDS
#ifdef WITH_IDS
          decodeTagID(rec->_tagid, id);

          sprintf(lines[l] + rec->_pos - winLo + 25, " %c "uint16FMTW(05)"-"uint16FMTW(05)"-"uint16FMTW(05)"-"uint16FMTW(05)" ",
                  (rec->_rev) ? '<' : '>',
                  id[0], id[1], id[2], id[3]);
#endif

          lineLen[l] = strlen(lines[l]);  //  Convert that trailing nul into a whitespace.
          lines[l][lineLen[l]] = ' ';

          uint32 err = 0;

          for (uint32 x=0; x<rec->_colorMismatch; x++) {
            uint32 pos = rec->_colorDiffs[err] & 0x3f;
            char   let = '*'; //bitsToColor[rec->_colorDiffs[err] >> 6];

            lines[l][rec->_pos - winLo + pos] = let;

            err++;
          }

          for (uint32 x=0; x<rec->_colorInconsistent; x++) {
            uint32 pos = rec->_colorDiffs[err] & 0x3f;
            char   let = bitsToColor[rec->_colorDiffs[err] >> 6];

            lines[l][rec->_pos - winLo + pos] = let;

            err++;
          }

          l = linesMax;
        }
        }
      }
    }

    bool stuff = false;

    for (uint32 i=0; i<linesMax; i++)
      if (lineLen[i] > 0)
        stuff = true;

    if (stuff) {
      fprintf(stdout, "\nALIGN "uint32FMT"-"uint32FMT"\n", winLo, winHi);
      fprintf(stdout, "     0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0\n");
      fprintf(stdout, "     012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890\n");

      for (uint32 i=0; i<linesMax; i++) {
        if (lineLen[i] > 0) {
          lines[i][lineLen[i]] = 0;
          fprintf(stdout, uint32FMTW(03)"] %s\n", i, lines[i]);
        }
      }
    }

    winLo = winHi;
    winHi = winLo + winSz;

    all.trimBeforeSeqPos(0, winLo);
  }

  delete inp;

  exit(0);
}











#if 0
    sprintf(linp, "rec "uint64HEX" "uint32FMT":"uint32FMT,
            rec->_tagid,
            rec->_seq,
            rec->_pos);
    while (*linp)
      linp++;

    uint32 err = 0;

    for (uint32 x=0; x<rec->_colorMismatch; x++) {
      sprintf(linp, " M:%c@%02d(%07d)",
              bitsToColor[rec->_colorDiffs[err] >> 6],
              (rec->_colorDiffs[err] & 0x3f),
              (rec->_colorDiffs[err] & 0x3f) + rec->_pos);
      while (*linp)
        linp++;
      err++;
    }

    for (uint32 x=0; x<rec->_colorInconsistent; x++) {
      sprintf(linp, " E:%c@%02d(%07d)",
              bitsToColor[rec->_colorDiffs[err] >> 6],
              (rec->_colorDiffs[err] & 0x3f),
              (rec->_colorDiffs[err] & 0x3f) + rec->_pos);
      while (*linp)
        linp++;
      err++;
    }

    fprintf(stdout, "%s\n", line);
#endif
