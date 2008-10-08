
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

#ifndef KMER_UTILS_H
#define KMER_UTILS_H

static const char *rcsid_KMER_UTILS_H = "$Id: kmer_utils.h,v 1.3 2008-10-08 22:02:57 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"

#ifndef minim
#define minim(a,b) ( (a) < (b) ? (a) : (b) )
#endif
#ifndef maxim
#define maxim(a,b) ( (a) > (b) ? (a) : (b) )
#endif


static int charval(char c){
  switch (c){
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
  default:
      return -1;
  }
}




static int kosherChar(char c){
  switch (c){
  case 'A':
  case 'a':
    break;
  case 'C':
  case 'c':
    break;
  case 'G':
  case 'g':
    break;
  case 'T':
  case 't':
    break;
  default:
    return 0;
  }
  return 1;
}


static uint64 oneBaseUpdateIdx(uint64 was, char nextch,int64 mask){
  uint64 rv=was;
  assert(kosherChar(nextch));
  rv <<= 2;
  rv &= mask;
  rv+= charval(nextch);
  return rv;
}


static void spellOutPacked(int idx, char *unpack,int unpackLen){
  int i,j;
  int mask= 3;
  assert(unpackLen<=sizeof(int)*4);
  unpack[unpackLen]='\0';
  for(i=0,j=unpackLen-1;i<unpackLen;i++,j--){
    int c = idx & mask;
    idx = (idx >> 2);
    switch (c){
    case 0:
      unpack[j] = 'A';
      break;
    case 1:
      unpack[j] = 'C';
      break;
    case 2:
      unpack[j] = 'G';
      break;
    case 3:
      unpack[j] = 'T';
      break;
    default:
      fprintf(stderr,"Unexpected condition packed char index = %d\n",c);
      exit(-1);
    }
  }
}





static void calc_kmer_profile(char *seq, int *kmerCounts, int merSize){
  int i;
  int N = strlen(seq);
  int kmers = (1<<(merSize*2));
  int32 goodchars=0;
  int64 merMask=0;
  uint32 idx=0;
  for(i=0;i<merSize;i++){
    merMask<<=2;
    merMask+=3;
  }


  for(i=0;i<kmers;i++){
    kmerCounts[i]=0;
  }

  for(i=0;i<N;i++){

    // check for unambig. DNA mer ... otherwise, continue
    if(!kosherChar(seq[i])){
      goodchars=0;
      idx=0;
      //      fprintf(stderr,"Skipping bad char at %d [%c]\n",i,seq[i]);
      continue;
    }
    goodchars++;
    if(goodchars > merSize)goodchars=merSize;
    idx = oneBaseUpdateIdx(idx,seq[i],merMask);
    if(goodchars!=merSize){
      //      fprintf(stderr,"Partial mer not fully processed, %d cur len = %d\n",
      //	      i,goodchars);
      continue;
    }
    kmerCounts[idx]++;
  }

}

static int rcKmer(int in,int ksize){
  int out=0;
  int i;

  for(i=0;i<ksize;i++){
    int bits = 3 - ( (in >> i*2 ) & 3);
    out += (bits << ((ksize-1-i)*2));
  }

  return out;

}
#endif
