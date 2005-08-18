#include <stdio.h>
#include "util++.H"
#include "constants.H"

#ifndef OVERLAP_H
#define OVERLAP_H

//  Define this to keep the sorted overlap file as ASCII -- it'll get
//  really big and be really slow, but you can then read it.
//
//#define ASCII_OVERLAPS

class  overlap_t {
public:
  u64bit   Aiid:26;  //  26 + 11 + 11 = 48
  u64bit   Abeg:11;
  u64bit   Aend:11;

  u64bit   ori:1;

  u64bit   Biid:26;
  u64bit   Bbeg:11;
  u64bit   Bend:11;

  u64bit   erate:16;

  void     decode(char *line, bool flip) {
    splitToWords  W(line);

    if (flip == false) {
      Aiid  = strtou32bit(W[0], 0L);
      Abeg  = strtou32bit(W[3], 0L);
      Aend  = strtou32bit(W[4], 0L);
      ori   = (W[2][0] == 'f');
      Biid  = strtou32bit(W[1], 0L);
      Bbeg  = strtou32bit(W[6], 0L);
      Bend  = strtou32bit(W[7], 0L);
      erate = (int)(atof(W[9]) * 100);
    } else {
      Aiid  = strtou32bit(W[1], 0L);
      Biid  = strtou32bit(W[0], 0L);
      erate = (int)(atof(W[9]) * 100);

      switch (W[2][0]) {
        case 'f':
          ori  = 1;
          Abeg = strtou32bit(W[6], 0L);
          Aend = strtou32bit(W[7], 0L);
          Bbeg = strtou32bit(W[3], 0L);
          Bend = strtou32bit(W[4], 0L);
          break;
        case 'r':
          ori  = 0;
          Abeg = strtou32bit(W[7], 0L);
          Aend = strtou32bit(W[6], 0L);
          Bbeg = strtou32bit(W[4], 0L);
          Bend = strtou32bit(W[3], 0L);
          break;
      } 
    }
  };


  bool     acceptable(void) {

    //  Remember if the 5' ends are far apart -- but we only care if
    //  it's a forward match.
    //
    bool far5prime = true;
    if (ori) {
      if (Abeg > Bbeg)
        far5prime = ((Abeg - Bbeg) > OBT_FAR5PRIME);
      else
        far5prime = ((Bbeg - Abeg) > OBT_FAR5PRIME);
    }

    //  Bbeg can greater than Bend, so we find the difference.  A is
    //  always properly ordered.
    //
    u32bit Adiff = Aend - Abeg;
    u32bit Bdiff = Bend - Bbeg;
    if (Bend < Bbeg)
      Bdiff = Bbeg - Bend;

    //  It's an acceptable overlap if the error is within tolerance, if
    //  it's long, and if the 5' ends are far apart.
    //
    //  We never use this anymore, the last version said erate <
    //  1.5, and both differences > 100.
    //
    return(((erate < OBT_MIN_ERATE) ||
            ((Adiff > OBT_MIN_DIFF) &&
             (Bdiff > OBT_MIN_DIFF))) &&
           (far5prime));
  };




  void     print(FILE *outf) const {
    u64bit  bogusLength = 666;
    fprintf(outf, u64bitFMTW(8)" "u64bitFMTW(8)" %c "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" %5.2f\n",
            Aiid, Biid,
            ori ? 'f' : 'r',
            Abeg, Aend, bogusLength,
            Bbeg, Bend, bogusLength,
            erate / 100.0);
  };

  void    dump(FILE *outf) const {
    fwrite(this, sizeof(overlap_t), 1, outf);
  };

  void    load(FILE *outf) {
    fread(this, sizeof(overlap_t), 1, outf);
  };

};

#if 0

inline
bool
acceptableOverlap(bool ori,
                  u32bit  leftA,
                  u32bit  rightA,
                  u32bit  leftB,
                  u32bit  rightB) {

  //  Remember if the 5' ends are far apart -- but we only care if
  //  it's a forward match.
  //
  bool far5prime = true;
  if (ori) {
    if (leftA > leftB)
      far5prime = ((leftA - leftB) > OBT_FAR5PRIME);
    else
      far5prime = ((leftB - leftA) > OBT_FAR5PRIME);
  }

  //  leftB can greater than Bend, so we find the difference.  A is
  //  always properly ordered.
  //
  u32bit Adiff = rightA - leftA;
  u32bit Bdiff = rightB - leftB;
  if (rightB < leftB)
    Bdiff = leftB - rightB;

  //  It's an acceptable overlap if the error is within tolerance, if
  //  it's long, and if the 5' ends are far apart.
  //
  //  We never use this anymore, the last version said erate <
  //  1.5, and both differences > 100.
  //
  bool  acceptable = (((erate < OBT_MIN_ERATE) ||
                       ((Adiff > OBT_MIN_DIFF) &&
                        (Bdiff > OBT_MIN_DIFF))) &&
                      (far5prime));

  return(acceptable);
}

#endif


#endif  //  OVERLAP_H
