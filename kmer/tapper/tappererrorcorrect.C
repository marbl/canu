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
    alignsLen      = new u32bit            [alignsMax];
    alignsPerBlock = 16384;
    alignsInp      = inp;

    for (u32bit i=0; i<alignsMax; i++) {
      aligns[i] = new tapperAlignment [alignsPerBlock];
      alignsLen[i] = alignsInp->getRecord(aligns[i], alignsPerBlock);
      fprintf(stderr, "block "u32bitFMT" has "u32bitFMT" things.\n", i, alignsLen[i]);
    }
  };

  ~alignmentList() {
    for (u32bit i=0; i<alignsMax; i++)
      delete [] aligns[i];

    delete [] aligns;
  };

  //  If the last element in the first block is below the specified
  //  seq,pos, we can dump all those alignments and get more.
  //
  void  trimBeforeSeqPos(u32bit seq, u32bit pos) {

  trimBeforeSeqPosAgain:
    if ((aligns[0][alignsLen[0]]._seq <= seq) &&
        (aligns[0][alignsLen[0]]._pos <  pos)) {
      tapperAlignment *save = aligns[0];

      for (u32bit i=1; i<alignsMax; i++) {
        aligns[i-1]    = aligns[i];
        alignsLen[i-1] = alignsLen[i];
      }

      aligns[alignsMax-1] = save;

      alignsLen[alignsMax-1] = alignsInp->getRecord(aligns[alignsMax-1], alignsPerBlock);

      fprintf(stderr, "block "u32bitFMT" has "u32bitFMT" things.\n", alignsMax-1, alignsLen[alignsMax-1]);

      goto trimBeforeSeqPosAgain;
    }
  };

  tapperAlignment *operator[](u32bit x) {
    u32bit block = x / alignsPerBlock;
    u32bit piece = x % alignsPerBlock;

    if (piece < alignsLen[block])
      return(aligns[block] + piece);

    return(0L);
  };

  bool empty(void) {
    return(alignsLen[0] == 0);
  };

private:
  u32bit            alignsMax;

  tapperAlignment **aligns;
  u32bit           *alignsLen;

  u32bit            alignsPerBlock;

  recordFile       *alignsInp;
};





int
main(int argc, char **argv) {
  char     *outputName  = 0L;
  char     *inputName   = 0L;

  u64bit    memoryLimit = 1024 * 1024 * 1024;

  {
  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-memory", 2) == 0) {
      memoryLimit = strtou64bit(argv[++arg], 0L) * 1024 * 1024;

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

  u32bit            winSz = 100;
  u32bit            winLo = 0;
  u32bit            winHi = winLo + winSz;

  char              lines[100][256];
  u32bit            lineLen[100];

  while (all.empty() == false) {
    memset(lines, ' ', sizeof(char) * 100 * 256);

    for (u32bit i=0; i<100; i++)
      lineLen[i] = 0;

    for (u32bit a=0; (all[a] != 0L) && (all[a]->_pos < winHi); a++) {
      tapperAlignment *rec = all[a];

      for (u32bit l=0; l<100; l++) {
        if ((winLo < rec->_pos) &&
            (lineLen[l] < rec->_pos - winLo)) {

          //fprintf(stdout, "at l="u32bitFMT" x="u32bitFMT"\n", l, rec->_pos - winLo);

          for (u32bit x=rec->_pos - winLo; x<rec->_pos - winLo + 25; x++)
            lines[l][x] = '.';

          lineLen[l] = rec->_pos - winLo + 25 + 1;

          u32bit err = 0;

          for (u32bit x=0; x<rec->_colorMismatch; x++) {
            u32bit pos = rec->_colorDiffs[err] & 0x3f;
            char   let = bitsToColor[rec->_colorDiffs[err] >> 6];

            lines[l][rec->_pos - winLo + pos] = let;

            err++;
          }

          for (u32bit x=0; x<rec->_colorInconsistent; x++) {
            u32bit pos = rec->_colorDiffs[err] & 0x3f;
            char   let = bitsToColor[rec->_colorDiffs[err] >> 6];

            lines[l][rec->_pos - winLo + pos] = let;

            err++;
          }

          l = 100;
        }
      }
    }

    fprintf(stdout, "\nALIGN "u32bitFMT"-"u32bitFMT"\n", winLo, winHi);

    for (u32bit i=0; i<100; i++) {
      if (lineLen[i] > 0) {
        lines[i][lineLen[i]  ] = 0;

        fprintf(stdout, u32bitFMTW(03)"] %s\n", i, lines[i]);
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
    sprintf(linp, "rec "u64bitHEX" "u32bitFMT":"u32bitFMT,
            rec->_tagid,
            rec->_seq,
            rec->_pos);
    while (*linp)
      linp++;

    u32bit err = 0;

    for (u32bit x=0; x<rec->_colorMismatch; x++) {
      sprintf(linp, " M:%c@%02d(%07d)",
              bitsToColor[rec->_colorDiffs[err] >> 6],
              (rec->_colorDiffs[err] & 0x3f),
              (rec->_colorDiffs[err] & 0x3f) + rec->_pos);
      while (*linp)
        linp++;
      err++;
    }

    for (u32bit x=0; x<rec->_colorInconsistent; x++) {
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
