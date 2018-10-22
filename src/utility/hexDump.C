
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
 *    Brian P. Walenz beginning on 2018-JUL-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"


//  Dump DATAlen bytes from DATA in a hex format.
//  It will print W bytes per line, separated into words of 8 bytes.
//  The end of the line will have the ASCII representation of the data.
//
//  00000000  00 01 02 03 04 05 06 07  08 09 0a 0b 0c 0d 0e 0f  '................'
//
void
hexDump(FILE *F,
        void *DATA, uint32 DATAlen, uint32 W) {
  char *STR = new char [8 + 2 + W * 3 + W/8 + 2 + 1 + W * 1 + 1];

  for (uint32 Dp=0; Dp < DATAlen; Dp += W) {
    uint8  *D   = (uint8 *)DATA + Dp;
    uint32  Ds  = (Dp + W <= DATAlen) ? (W) : (DATAlen - Dp);

    for (uint32 Z=Dp, ii=8; ii>0; Z>>=4)  //  Dump the address in hexadecimal
      STR[--ii] = ((Z & 0x0f) < 0x0a) ? ((Z & 0x0f) + '0') : ((Z & 0x0f) - 0x0a + 'a');

    char *H = STR + 8;                     //  Data pointer
    char *A = STR + 8 + 1 + 3 * W + W/8;   //  ASCII pointer

    *H++  = ' ';                      //  Another space is added at ii=0 below.

    *A++  = ' ';                      //  Space between the last digit and the string
    *A++ = '\'';                      //  Bracket at the start of the string.

    for (uint32 ii=0; ii<W; ii++) {
      if ((ii % 8) == 0)              //  An extra space between words
        *H++ = ' ';

      if (ii < Ds) {                  //  Emit a digit, or...
        *H++ = ((D[ii] & 0xf0) < 0xa0) ? (((D[ii] & 0xf0) >> 4) + '0') : (((D[ii] & 0xf0) >> 4) - 0x0a + 'a');
        *H++ = ((D[ii] & 0x0f) < 0x0a) ? (((D[ii] & 0x0f)     ) + '0') : (((D[ii] & 0x0f)     ) - 0x0a + 'a');
      }
      else {
        *H++ = ' ';                   //  ...spaces if we fell off the end of the data
        *H++ = ' ';
      }

      *H++ = ' ';                     //  Space between digits

      if (ii < Ds)                    //  Printable ASCII or a dot
        *A++ = ((' ' <= D[ii]) && (D[ii] <= '~')) ? (D[ii]) : ('.');
    }

    *A++ = '\'';                      //  Bracket at the end of the string.
    *A++ = '\n';

    *A = 0;                           //  NUL terminate the string.

    fputs(STR, F);
  }

  delete [] STR;
}


