#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util++.H"
#include "constants.H"

#ifndef OVERLAP_H
#define OVERLAP_H


//  Define this if you are going to be cross-platform or
//  cross-architecture -- bitfields are not guarateed to be binary
//  compatible.
#define SAFE


class  overlap_t {
public:
  u64bit   Aiid:26;
  u64bit   Abeg:11;
  u64bit   Aend:11;
  u64bit   ori:1;
  u64bit   pad:15;

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
      erate = (int)(0.5 + atof(W[9]) * 100.0);  //  0.5 to prevent stupid boundary cases in truncation
    } else {
      Aiid  = strtou32bit(W[1], 0L);
      Biid  = strtou32bit(W[0], 0L);
      erate = (int)(0.5 + atof(W[9]) * 100.0);  //  0.5 to prevent stupid boundary cases in truncation

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


  void     flip(void) {
    u64bit t;
    t = Aiid; Aiid = Biid; Biid = t;
    if (ori) {
      //  Forward, so just swap the A and B coords
      t=Abeg; Abeg=Bbeg; Bbeg=t;
      t=Aend; Aend=Bend; Bend=t;
    } else {
      //  Reverse, swap the A and B coords, then reverse both
      t=Abeg; Abeg=Bend; Bend=t;
      t=Aend; Aend=Bbeg; Bbeg=t;
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
    //  We never use this anymore: the last version said erate <
    //  1.5, and both differences > 100.
    //
    return(((erate < OBT_MIN_ERATE) ||
            ((Adiff > OBT_MIN_DIFF) &&
             (Bdiff > OBT_MIN_DIFF))) &&
           (far5prime));
  };




  void     print(FILE *outf) const {
    u64bit  bogusLength = 666;
    fprintf(outf, u64bitFMTW(8)" "u64bitFMTW(8)" %c "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(5)"."u64bitFMTW(02)"\n",
            Aiid, Biid,
            ori ? 'f' : 'r',
            Abeg, Aend, bogusLength,
            Bbeg, Bend, bogusLength,
            erate / 100, erate % 100);
  };

  void    dump(FILE *outf, bool fast=false) const {

#ifndef SAFE
    fast = true;
#endif

    if (fast) {
      fwrite(this, sizeof(overlap_t), 1, outf);
    } else {
      u64bit a[2] = {u64bitZERO, u64bitZERO};
      a[0]  = Aiid;     a[0] <<= 11;
      a[0] |= Abeg;     a[0] <<= 11;
      a[0] |= Aend;     a[0] <<= 1;
      a[0] |= ori;

      a[1]  = Biid;     a[1] <<= 11;
      a[1] |= Bbeg;     a[1] <<= 11;
      a[1] |= Bend;     a[1] <<= 16;
      a[1] |= erate;
      fwrite(a, sizeof(u64bit), 2, outf);
    }
  };

  void    load(FILE *outf, bool fast=false) {

#ifndef SAFE
    fast = true;
#endif

    if (fast) {
      fread(this, sizeof(overlap_t), 1, outf);
    } else {
      u64bit a[2];
      fread(a, sizeof(u64bit), 2, outf);
      ori   = a[0] & u64bitMASK(1);   a[0] >>= 1;
      Aend  = a[0] & u64bitMASK(11);  a[0] >>= 11;
      Abeg  = a[0] & u64bitMASK(11);  a[0] >>= 11;
      Aiid  = a[0] & u64bitMASK(26);

      erate = a[1] & u64bitMASK(16);  a[1] >>= 16;
      Bend  = a[1] & u64bitMASK(11);  a[1] >>= 11;
      Bbeg  = a[1] & u64bitMASK(11);  a[1] >>= 11;
      Biid  = a[1] & u64bitMASK(26);
    }
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
