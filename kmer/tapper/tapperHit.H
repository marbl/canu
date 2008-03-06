class tapperHit {
public:
  //tapperHit() {};
  //~tapperHit();

  void writeHit(FILE *OF,
                char *tag1name, kMer& tag1f, kMer& tag1r,
                char *tag2name, kMer& tag2f, kMer& tag2r) {

    //  Hits come out as (all one line, tab separated):
    //
    //  refSeqId
    //  tag1Id tag1pos tag1/ref
    //  ori
    //  tag1Id tag1pos tag1/ref

    char  tag1seq[256] = {0};
    char  ref1seq[256] = {0};

    char  tag2seq[256] = {0};
    char  ref2seq[256] = {0};

    strcpy(tag1seq, "1nothing1");
    strcpy(ref1seq, "1NOTHING1");
    strcpy(tag2seq, "2nothing2");
    strcpy(ref2seq, "2NOTHING2");

    char  ori1[4] = {'-', '-', 0, 0};
    char  ori2[4] = {'-', '-', 0, 0};

    if ((ori == 0) || (ori == 4) || (ori == 5))
      ori1[1] = '>';
    if ((ori == 1) || (ori == 6) || (ori == 7))
      ori1[0] = '<';

    if ((ori == 2) || (ori == 4) || (ori == 6))
      ori2[1] = '>';
    if ((ori == 3) || (ori == 5) || (ori == 7))
      ori2[0] = '<';

    switch (ori) {
      case 0:
      case 1:
        fprintf(OF, u64bitFMT"\t%s\t"u64bitFMT"\t%s/%s\t%s\n",
                seqIdx,
                tag1name, tag1pos, tag1seq, ref1seq, ori1);
        break;
      case 2:
      case 3:
        fprintf(OF, u64bitFMT"\t%s\t"u64bitFMT"\t%s/%s\t%s\n",
                seqIdx,
                tag2name, tag2pos, tag2seq, ref2seq, ori2);
        break;
      case 4:
      case 5:
      case 6:
      case 7:
        fprintf(OF, u64bitFMT"\t%s\t"u64bitFMT"\t%s/%s\t%s\t%s\t"u64bitFMT"\t%s/%s\t%s\n",
                seqIdx,
                tag1name, tag1pos, tag1seq, ref1seq, ori1,
                tag2name, tag2pos, tag2seq, ref2seq, ori2);
        break;
    }
  };

  //  ori  tag1  tag2  t1m t1r t2m t2r
  //
  //    0   -->   x    10 0* --  8 &  9
  //    1   <--   x    11 0* -- 12 & 13
  //    2    x   -->   0* 10 --  2 & 6
  //    3    x   <--   0* 11 --  3 & 7
  //    4   -->  -->   10 10 -- 10
  //    5   -->  <--   10 11 -- 11
  //    6   <--  -->   11 10 -- 14
  //    7   <--  <--   11 11 -- 15
  //
  //  Invalid values -- 0, 1, 4, 5
  //
  void   setOrientation(u32bit tag1mapped, u32bit tag1reversed,
                        u32bit tag2mapped, u32bit tag2reversed) {

    //                      0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15
    u32bit  oLookup[16] = { 0, 0, 2, 3, 0, 0, 2, 3, 0, 0, 4, 5, 1, 1, 6, 7};
    ori = oLookup[(tag1mapped << 3) | (tag1reversed << 2) | (tag2mapped << 1) | (tag2reversed)];
  };


#warning need a real consturctor that sets the errors found.

public:
  //  Output is a list of positions we map to.  255 bits -> 32 bytes per tag hit.

  u64bit  seqIdx      : 24;  //  index of the reference sequence (up to 16 million)
  u64bit  tag1errors  : 4;   //  tag 1 errors are stored first,
  u64bit  tag2errors  : 4;   //  then tag 2 errors are stored
  u64bit  changes     : 32;  //  16 errors, 2 bits to tell reference ACGTN

  u64bit  unused      :  1;
  u64bit  ori         :  3;  //  4 possible single tag orientations, 4 possible two tag orientations
  u64bit  tag1pos     : 30;  //  references up to 1gbp
  u64bit  tag2pos     : 30;  //  references up to 1gbp

  u8bit   positions[16];     //  8 bits - position the error occurs at in the tag
};

