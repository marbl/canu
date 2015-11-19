
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
 *    Brian P. Walenz from 2015-FEB-09 to 2015-JUL-08
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "Display_Alignment.H"
#include "gkStore.H"

#define  DISPLAY_WIDTH   250

//  Show (to  stdout ) the alignment encoded in  delta [0 .. (delta_ct - 1)]
//  between strings  a [0 .. (a_len - 1)]  and  b [0 .. (b_len - 1)] .

void
Display_Alignment(char  *a,    int32   a_len,
                  char  *b,    int32   b_len,
                  int32 *delta,
                  int32  delta_ct) {

  char  *top = new char [AS_MAX_READLEN + 1];
  char  *bot = new char [AS_MAX_READLEN + 1];

  int32 nAgap   = 0;
  int32 top_len = 0;

  {
    int32 i = 0;
    int32 j = 0;

    for (int32 k = 0;  k < delta_ct;  k++) {
      for (int32 m = 1;  m < abs (delta[k]);  m++) {
        top[top_len++] = a[i++];
        j++;
      }

      if (delta[k] < 0) {
        top[top_len++] = '-';
        j++;
        nAgap++;

      } else {
        top[top_len++] = a[i++];
      }
    }

    while (i < a_len && j < b_len) {
      top[top_len++] = a[i++];
      j++;
    }
    top[top_len] = '\0';
  }


  int32 nBgap   = 0;
  int32 bot_len = 0;

  {
    int32 i = 0;
    int32 j = 0;

    for (int32 k = 0;  k < delta_ct;  k++) {
      for (int32 m = 1;  m < abs (delta[k]);  m++) {
        bot[bot_len++] = b[j++];
        i++;
      }

      if (delta[k] > 0) {
        bot[bot_len++] = '-';
        i++;
        nBgap++;

      } else {
        bot[bot_len++] = b[j++];
      }
    }

    while (j < b_len && i < a_len) {
      bot[bot_len++] = b[j++];
      i++;
    }

    bot[bot_len] = '\0';
  }


  uint32 diffs = 0;

  for (int32 i=0; (i < top_len) || (i < bot_len); i += DISPLAY_WIDTH) {
    fprintf(stderr, "%d\n", i);
    fprintf(stderr, "A: ");

    for (int32 j=0;  (j < DISPLAY_WIDTH) && (i+j < top_len);  j++)
      putc(top[i+j], stderr);

    fprintf(stderr, "\n");
    fprintf(stderr, "B: ");

    for (int32 j=0; (j < DISPLAY_WIDTH) && (i+j < bot_len); j++)
      putc(bot[i+j], stderr);

    fprintf(stderr, "\n");
    fprintf(stderr, "   ");

    for (int32 j=0; (j<DISPLAY_WIDTH) && (i+j < bot_len) && (i+j < top_len); j++)
      if ((top[i+j] != ' ') && (bot[i+j] != ' ') && (tolower(top[i+j]) != tolower(bot[i+j]))) {
        diffs++;
        putc('^', stderr);
      } else {
        putc(' ', stderr);
      }

    putc('\n', stderr);
  }

  fprintf(stderr, "differences: %u   Agaps: %d   Bgaps: %d)\n", diffs, nAgap, nBgap);

  delete [] top;
  delete [] bot;
}
