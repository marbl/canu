
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
 *  This file is derived from:
 *
 *    src/stores/gkStoreEncode.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2015-DEC-03
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sqStore.H"


//  Encode seq as 2-bit bases.  Doesn't touch qlt.
uint32
sqReadData::sqReadData_encode2bit(uint8 *&chunk, char *seq, uint32 seqLen) {

  //  Scan the read, if there are non-acgt, return length 0; this cannot encode it.

  for (uint32 ii=0; ii<seqLen; ii++) {
    char  base = seq[ii];

    if ((base != 'a') && (base != 'A') &&
        (base != 'c') && (base != 'C') &&
        (base != 'g') && (base != 'G') &&
        (base != 't') && (base != 'T'))
      return(0);
  }

  uint8  acgt[256] = { 0 };

  acgt['a'] = acgt['A'] = 0x00;
  acgt['c'] = acgt['C'] = 0x01;
  acgt['g'] = acgt['G'] = 0x02;
  acgt['t'] = acgt['T'] = 0x03;

  uint32 chunkLen = 0;

  if (chunk == NULL)
    chunk = new uint8 [ seqLen / 4 + 1];

  for (uint32 ii=0; ii<seqLen; ) {
    uint8  byte = 0;

    if (ii + 4 < seqLen) {
      byte  = acgt[seq[ii++]];  byte <<= 2;
      byte |= acgt[seq[ii++]];  byte <<= 2;
      byte |= acgt[seq[ii++]];  byte <<= 2;
      byte |= acgt[seq[ii++]];
    }

    else {
      if (ii < seqLen)  { byte |= acgt[seq[ii++]]; }   byte <<= 2;  //  The if here is redundant, but pretty.
      if (ii < seqLen)  { byte |= acgt[seq[ii++]]; }   byte <<= 2;  //  Yes, everything shifts, not a mistake to leave out the braces.
      if (ii < seqLen)  { byte |= acgt[seq[ii++]]; }   byte <<= 2;
      if (ii < seqLen)  { byte |= acgt[seq[ii++]]; }
    }

    chunk[chunkLen++] = byte;
  }

  return(chunkLen);
}



bool
sqReadData::sqReadData_decode2bit(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {

  if (chunkLen == 0)
    return(false);

  uint32   chunkPos = 0;

  char     acgt[4] = { 'A', 'C', 'G', 'T' };

  for (uint32 ii=0; ii<seqLen; ) {
    assert(chunkPos < chunkLen);

    uint8  byte = chunk[chunkPos++];

    if (ii + 4 < seqLen) {
      seq[ii++] = acgt[((byte >> 6) & 0x03)];
      seq[ii++] = acgt[((byte >> 4) & 0x03)];
      seq[ii++] = acgt[((byte >> 2) & 0x03)];
      seq[ii++] = acgt[((byte >> 0) & 0x03)];
    }

    else {
      if (ii < seqLen)  seq[ii++] = acgt[((byte >> 6) & 0x03)];  //  This if is also redundant, and also pretty.
      if (ii < seqLen)  seq[ii++] = acgt[((byte >> 4) & 0x03)];
      if (ii < seqLen)  seq[ii++] = acgt[((byte >> 2) & 0x03)];
      if (ii < seqLen)  seq[ii++] = acgt[((byte >> 0) & 0x03)];
    }
  }

  seq[seqLen] = 0;

  return(true);
}



//  Encode seq as 3-bases-in-7-bits.  Doesn't touch qlt.
uint32
sqReadData::sqReadData_encode3bit(uint8 *&UNUSED(chunk), char *UNUSED(seq), uint32 UNUSED(seqLen)) {
  return(0);
}

bool
sqReadData::sqReadData_decode3bit(uint8 *UNUSED(chunk), uint32 UNUSED(chunkLen), char *UNUSED(seq), uint32 UNUSED(seqLen)) {
  return(false);
}





//  Encode qualities as 4 bit integers.  Doesn't touch seq.
uint32
sqReadData::sqReadData_encode4bit(uint8 *&UNUSED(chunk), uint8 *UNUSED(qlt), uint32 UNUSED(qltLen)) {
  return(0);
}

bool
sqReadData::sqReadData_decode4bit(uint8 *UNUSED(chunk), uint32 UNUSED(chunkLen), uint8 *UNUSED(qlt), uint32 UNUSED(qltLen)) {
  return(false);
}





//  Encode qualities as 5 bit integers.  Doesn't touch seq.
uint32
sqReadData::sqReadData_encode5bit(uint8 *&UNUSED(chunk), uint8 *UNUSED(qlt), uint32 UNUSED(qltLen)) {
  return(0);
}

bool
sqReadData::sqReadData_decode5bit(uint8 *UNUSED(chunk), uint32 UNUSED(chunkLen), uint8 *UNUSED(qlt), uint32 UNUSED(qltLen)) {
  return(false);
}



//  Analyze a QV string to decide if all the QVs are the same.  Returns a valid QV if so, 255 otherwise.
//
//  This occurs when we store a default QV, then load the read (which expands the default QV to the
//  QV array) and we then store the sequence again.  As when we load corrected reads.
//
uint32
sqReadData::sqReadData_encodeConstantQV(uint8 *qlt, uint32 qltLen, uint32 defaultQV) {

  if (qltLen == 0)                     //  If no sequence, do nothing.
    return(255);

  if (qlt[0] == 255)                   //  If no QV's, return the
    return(defaultQV);                 //  default QV for the library.

  for (uint32 ii=1; ii<qltLen; ii++)
    if (qlt[0] != qlt[ii])             //  If the QV's differ,
      return(255);                     //  return the sentinel.

  return(qlt[0]);                      //  Otherwise, QV's are the same!
}
