 
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sequence/sequence.H"



void
doShiftRegister(shiftRegisterParameters &srPar) {
  uint32 len = 0;
  char   sr[17];
  uint8  sv[17];

  srPar.initialize();

  len = srPar.len;
  len = 8;

  for (uint32 ii=0; ii<srPar.len; ii++) {
    sr[ii] =  srPar.sr[ii];
    sv[ii] = (srPar.sv[ii] == '1') ? 1 : 0;
  }

  sr[0] = 'A';  sv[0] = 1;   //  Oldest value
  sr[1] = 'A';  sv[1] = 0;
  sr[2] = 'A';  sv[2] = 0;
  sr[3] = 'A';  sv[3] = 0;
  sr[4] = 'A';  sv[4] = 0;
  sr[5] = 'A';  sv[5] = 0;
  sr[6] = 'A';  sv[6] = 0;
  sr[7] = 'G';  sv[7] = 0;  //  Newest value

  sr[srPar.len] = 0;
  sv[srPar.len] = 0;

  uint32   kmer   = 0;
  uint32  *detect = new uint32 [65536];

  while ((sv[7] != 0) ||
         (sv[6] != 0) ||
         (sv[5] != 0) ||
         (sv[4] != 0) ||
         (sv[3] != 0) ||
         (sv[2] != 0) ||
         (sv[1] != 0) ||
         (sv[0] != 0)) {

    sr[0] = 'A';
    sr[1] = 'A';
    sr[2] = 'A';
    sr[3] = 'A';
    sr[4] = 'A';
    sr[5] = 'A';
    sr[6] = 'A';
    sr[7] = 'G';

    kmer = 0x0003;

#define STRING
#undef  KMER

#ifdef STRING
      fprintf(stdout, "%s", sr);
#endif

    memset(detect, 0, sizeof(uint32) * 65536);
    detect[kmer] = 1;

    for (uint32 ii=0; ii<65536; ii++) {
#ifdef KMER
      fprintf(stdout, "%06u %s 0x%04x %u\n", ii, sr, kmer, detect[kmer]);
#endif

      uint32  next = 0;

      for (uint32 kk=0; kk<len; kk++)
        if (sv[kk])
          next += (sr[kk] >> 1) & 0x03;

      if      ((next & 0x03) == 0x00)  next = 'A';
      else if ((next & 0x03) == 0x01)  next = 'C';
      else if ((next & 0x03) == 0x03)  next = 'G';
      else if ((next & 0x03) == 0x02)  next = 'T';

      for (uint32 kk=0; kk<len-1; kk++)
        sr[kk] = sr[kk+1];

      sr[len-1] = next;

#ifdef STRING
      fprintf(stdout, "%c", next);
#endif

      kmer <<= 2;
      kmer  |= ((sr[len-1] >> 1) & 0x3);
      kmer &= 0xffff;

      detect[kmer]++;

      if (detect[kmer] == 2) {
        if (ii > 0) {
#ifdef KMER
          fprintf(stdout, "%06u %s 0x%04x %u\n", ii, sr, kmer, detect[kmer]);
#endif
#ifdef STRING
          fprintf(stdout, "\n");
#endif
          fprintf(stderr, "cycle at ii=%5u  %d%d%d%d%d%d%d%d\n",
                  ii, sv[0], sv[1], sv[2], sv[3], sv[4], sv[5], sv[6], sv[7]);
        }
        break;
      }
    }

    sv[7]++;

    if (sv[7] > 1)  { sv[7] = 0;  sv[6]++; }
    if (sv[6] > 1)  { sv[6] = 0;  sv[5]++; }
    if (sv[5] > 1)  { sv[5] = 0;  sv[4]++; }
    if (sv[4] > 1)  { sv[4] = 0;  sv[3]++; }
    if (sv[3] > 1)  { sv[3] = 0;  sv[2]++; }
    if (sv[2] > 1)  { sv[2] = 0;  sv[1]++; }
    if (sv[1] > 1)  { sv[1] = 0;  sv[0]++; }
    if (sv[0] > 1)  { break; }
  }

  delete [] detect;
}


