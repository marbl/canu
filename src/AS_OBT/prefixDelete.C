
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

#include "util++.H"
#include "trim.H"

//  Removes all reads that are a perfect prefix of some other read.
//  Reads must be in a library that has the deletePerfectPrefixes flag
//  set.
//
//  The algorithm builds a 64-bit value from the first N bases (N <=
//  AS_FRAG_MIN_LEN, currently 48 <= 64), sorts the hashes, then
//  examines any clique of hash collisions for perfect prefixes.


struct fragHash {
  uint64    hash;
  uint32    iid;
};


int
fragHashCompare(const void *a, const void *b) {
  fragHash const *A = (fragHash const *)a;
  fragHash const *B = (fragHash const *)b;

  if (A->hash < B->hash) return(-1);
  if (A->hash > B->hash) return( 1);
  return(0);
}


int
main(int argc, char **argv) {
  char             *gkpName = 0L;
  FILE             *logFile = 0L;
  GateKeeperStore  *gkp     = 0L;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-frg", 2) == 0) {
      gkpName = argv[++arg];
    } else if (strncmp(argv[arg], "-log", 2) == 0) {
      errno=0;
      logFile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open %s for writing the log: %s\n", argv[arg], strerror(errno)), exit(1);
    } else {
      fprintf(stderr, "Invalid option: '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err) || (gkpName == 0L)) {
    fprintf(stderr, "usage: %s -frg some.gkpStore\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -frg       Operate on this gkpStore\n");
    fprintf(stderr, "  -log       Report the iid, original trim and new quality trim\n");
    exit(1);
  }

  gkp = openGateKeeperStore(gkpName, TRUE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpName);
    exit(1);
  }

  gkp->frg = convertStoreToMemoryStore(gkp->frg);

  uint32        firstElem = getFirstElemFragStore(gkp);
  uint32        lastElem  = getLastElemFragStore(gkp) + 1;

  fragRecord    fr1;
  fragRecord    fr2;

  fragHash   *fh    = new fragHash [lastElem - firstElem + 1];
  uint32      fhLen = 0;

  fprintf(stderr, "Store contains "F_U32" fragments.\n", lastElem - firstElem + 1);

  uint64 map[256] = { 0 };
  uint32 s, h, n;

  for (s=0; s<256; s++)
    map[s] = 0;
  map['A'] = map['a'] = 0x00;
  map['C'] = map['c'] = 0x01;
  map['G'] = map['g'] = 0x02;
  map['T'] = map['t'] = 0x03;

  for (uint32 elem=firstElem; elem<lastElem; elem++) {

    getFrag(gkp, elem, &fr1, FRAG_S_INF);

    if (getFragRecordIsDeleted(&fr1))
      continue;

    //  Check if we can delete prefixes from this read.  Assumes the
    //  gkpStore caches libraries, which it does (or did).
    //
    {
      GateKeeperLibraryRecord  *gkpl = getGateKeeperLibrary(gkp, getFragRecordLibraryIID(&fr1));
      if ((gkpl == NULL) || (gkpl->deletePerfectPrefixes == 0))
        continue;
    }

    getFrag(gkp, elem, &fr1, FRAG_S_SEQ);

    char *seq1 = getFragRecordSequence(&fr1);

    uint32 seqLen   = getFragRecordSequenceLength(&fr1);
    uint64 hash     = 0;

    assert(seqLen >= 48);

    //  Our "hash" is just the spaced seed "101" (repeating).  It
    //  covers the first 48 bases, picking out 32.
    //
    for (s=0, n=0; n<16; n++) {
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
      s++;
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
    }

    fh[fhLen].hash     = hash;
    fh[fhLen].iid      = elem;

    fhLen++;
  }

  fprintf(stderr, "Store contains "F_U32" fragments allowing exact prefixes to be deleted.\n", fhLen);

  qsort(fh, fhLen, sizeof(fragHash), fragHashCompare);

  uint32  beg = 0;
  uint32  end = 0;

  while (beg < fhLen) {

    //  We DO need to examine the whole clique (pairwise).  We cannot
    //  simply sort by size, because if we get three frags of the same
    //  size, it could be that #1 is a prefix of #3, and #2 is just of
    //  the same size.  Even there, we'd need to examine all pairs.

    end = beg + 1;

    //  First, find a pair of adjacent matches
    //
    while ((end < fhLen) &&
           (fh[beg].hash != fh[end].hash)) {
      beg++;
      end++;
    }

    //  Got a match?
    //
    if (end < fhLen) {

      //  Advance end to the end of the matches
      //
      while ((fh[beg].hash == fh[end].hash) && (end < fhLen))
        end++;

      //  Yeah, we could extend scope of this test to include the for
      //  loops, but those will stop quick enough.

      //if (beg + 1 < end)
      //  fprintf(stderr, "Clique from "F_U32" to "F_U32" ("F_U32" things)\n", beg, end, end - beg);

      //  Compare all-vs-all in the range
      //
      for (uint32 b=beg; b<end; b++) {
        for (uint32 e=b+1; e<end; e++) {

          AS_IID     iid1 = fh[b].iid;
          AS_IID     iid2 = fh[e].iid;

          getFrag(gkp, iid1, &fr1, FRAG_S_SEQ);
          getFrag(gkp, iid2, &fr2, FRAG_S_SEQ);

          uint32 del1 = getFragRecordIsDeleted(&fr1);
          uint32 del2 = getFragRecordIsDeleted(&fr2);

          uint32 len1 = getFragRecordSequenceLength(&fr1);
          uint32 len2 = getFragRecordSequenceLength(&fr2);

          if ((del1) && (len1 < len2))
            continue;
          if ((del2) && (len2 < len1))
            continue;

          if (len1 == len2) {
            if ((del1) && (iid1 < iid2))
              continue;
            if ((del2) && (iid2 < iid1))
              continue;
          }

          if (del1 && del2)
            continue;

          char *seq1 = getFragRecordSequence(&fr1);
          char *seq2 = getFragRecordSequence(&fr2);

          uint32 len = MIN(len1, len2);

          if (strncmp(seq1, seq2, len) == 0) {

            //  A real collision.  Delete smaller of the two (either
            //  smaller sequence length or smaller iid).  We can skip
            //  the delete if it's already deleted.

            AS_UID     deletedUID = AS_UID_undefined();
            AS_IID     deletedIID = 0;
            uint32     deleted    = 0;

            if ((len == getFragRecordSequenceLength(&fr1)) &&
                (len == getFragRecordSequenceLength(&fr2))) {
              deletedIID = (iid1 < iid2) ? iid1 : iid2;
              deletedUID = (iid1 < iid2) ? getFragRecordUID(&fr1) : getFragRecordUID(&fr2);
              deleted    = (iid1 < iid2) ? del1 : del2;
            } else if (len == getFragRecordSequenceLength(&fr1)) {
              deletedIID = iid1;
              deletedUID = getFragRecordUID(&fr1);
              deleted    = del1;
            } else {
              deletedIID = iid2;
              deletedUID = getFragRecordUID(&fr2);
              deleted    = del2;
            }

            if (deleted == 0) {
              delFrag(gkp, deletedIID);
              fprintf(stdout, "%s,%d\t%s,%d\tlen="F_U32"\tdeleted=%s\n",
                      AS_UID_toString1(getFragRecordUID(&fr1)), getFragRecordIID(&fr1),
                      AS_UID_toString2(getFragRecordUID(&fr2)), getFragRecordIID(&fr2),
                      len,
                      AS_UID_toString(deletedUID));
            }
          }
        }
      }
    }

    beg = end;
  }

  closeGateKeeperStore(gkp);
}


