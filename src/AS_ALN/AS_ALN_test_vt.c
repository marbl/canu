
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
/*********************************************************************/
// headers
/*********************************************************************/
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>

// project headers
#include "AS_global.h"
#include "AS_ALN_aligners.h"

char Bases[4] = {'A','C','G','T'};

cds_int8 BaseComplement[256] =
{
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
/* @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
   P   Q   R   S   T   U   V   W   X   Y   Z                */
  'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
  'N','N','N','N','A','A','N','N','N','N','N','N','N','N','N','N',
  'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
  'N','N','N','N','A','A','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
};

// ASCII to compressed bit representation for nucleotides
cds_int8 BitEquivalent[256] =
{
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
/* @  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O
   P  Q  R  S  T  U  V  W  X  Y  Z                */
  -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1, 3, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1, 3, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
};

static void MuckUpSequence( char * seq, cds_float32 var )
{
  int   i, j;
  int   num_subindels = 0;
  float prob;
  int   subindel;
  
  for( i = 0; seq[i] != '\0'; i++ )
  {
    if( (prob = ((float) rand()) / RAND_MAX) < var )
    {
      num_subindels++;
      subindel = rand() % 3;
      switch( subindel )
      {
        case 0:
          // substitution
          seq[i] = Bases[(BitEquivalent[(int) seq[i]]+1)%4];
          break;
        case 1:
          // insertion
          for( j = strlen( seq ) + 1; j > i; j-- )
            seq[j] = seq[j-1];
          seq[i++] = Bases[random() % 4];
          break;
        case 2:
          // deletion
          for( j = i; seq[j] != '\0'; j++ )
            seq[j] = seq[j+1];
          i--;
          break;
      }
    }
  }
  //  fprintf( stderr, "var: %f, pct indels: %d\n", var, num_subindels * 100 / strlen( seq ) );
}


static void ReverseComplement( char * dest, char * src )
{
  int i;
  int len = strlen( src );
  for( i = 0; i < len; i++ )
    dest[i] = BaseComplement[(int)src[len - 1 - i]];
  dest[i] = '\0';
}


typedef struct
{
  SeqInterval   seq_interval;
} TestSequence;

typedef struct
{
  SeqInterval   frag;
  SeqInterval   frag_clr;
  DirectionType frag_dir;
  SeqInterval   vector;
  VectorType    vector_type;
  SeqInterval   answer;
  cds_float32   variation;
  cds_int32     min_length;
} TestCase;

#define NUM_TEST_CASES  8
TestCase Cases[NUM_TEST_CASES] =
{
  /*       10        490
    frag: 0------------>500
    vect:         460---------->1000
     ans:  10-----460
  */
  {{   0, 500},{  10, 490},AS_FORWARD,{ 460,1000},AS_INSERT_VECTOR,{  10, 460},.10,10},
  {{ 600,1000},{  10, 350},AS_REVERSE,{ 100, 680},AS_INSERT_VECTOR,{  10, 320},.10,10},
  {{ 600,1000},{  10, 350},AS_REVERSE,{ 100, 670},AS_INSERT_VECTOR,{  10, 330},.10,10},
  {{ 600,1000},{  10, 350},AS_REVERSE,{ 100, 660},AS_INSERT_VECTOR,{  10, 340},.10,5},
  {{ 600,1000},{  10, 350},AS_REVERSE,{ 100, 654},AS_INSERT_VECTOR,{  10, 346},.10,3},
  {{ 600,1000},{  10, 350},AS_REVERSE,{ 100, 640},AS_INSERT_VECTOR,{  10, 350},.10,3},
  {{ 600,1000},{  10, 390},AS_FORWARD,{ 100, 650},AS_INSERT_VECTOR,{  10, 390},.10,3},
  {{ 600,1000},{  10, 390},AS_REVERSE,{ 900,1100},AS_INSERT_VECTOR,{  10, 390},.10,3}
};


static void CutSequence( char * seq, char * bases, SeqInterval interval, DirectionType dir )
{
  char   rev[AS_READ_MAX_LEN];
  char * temp_seq;

  temp_seq = (dir == AS_REVERSE) ? rev : seq;
  
  strncpy( temp_seq, &(bases[interval.bgn]), interval.end - interval.bgn );
  temp_seq[interval.end - interval.bgn] = '\0';
  
  if( dir == AS_REVERSE )
    ReverseComplement( seq, temp_seq );
}


int main( int argc, char ** argv )
{
  InternalFragMesg ifg;
  char frag_sequence[AS_READ_MAX_LEN];
  InternalScreenItemMesg isn;
  char vector_sequence[AS_READ_MAX_LEN];
  char bases[AS_READ_MAX_LEN];
  int i, k;

  // create a random sequence to seed frags & vectors
  srandom( time( 0 ) );
  for( i = 0; i < AS_READ_MAX_LEN - 1; i++ )
    bases[i] = Bases[random() % 4];
  bases[i] = '\0';

  // set unchanging elements of ifg & isn
  ifg.eaccession = ifg.iaccession = 1;
  ifg.action = AS_ADD;
  ifg.type = AS_READ;
  ifg.sequence = frag_sequence;
  ifg.quality = NULL;
  ifg.source = "\0";
  ifg.screened = NULL;
  
  isn.eaccession = isn.iaccession = 2;
  isn.action = AS_ADD;
  isn.sequence = vector_sequence;
  isn.source = "\0";

  // loop through testcases
  for( k = 0; k < 100; k++ )
  {
  for( i = 0; i < NUM_TEST_CASES; i++ )
  {
    // cut out a fragment sequence
    CutSequence( ifg.sequence, bases, Cases[i].frag, Cases[i].frag_dir );
    MuckUpSequence( ifg.sequence, 0.02 /*Cases[i].variation*/ );
    ifg.clear_rng.bgn = Cases[i].frag_clr.bgn;
    ifg.clear_rng.end = Cases[i].frag_clr.end;

    // cut out a vector sequence
    CutSequence( isn.sequence, bases, Cases[i].vector, AS_FORWARD );

    // Trim the fragment
    DP_Trim_Vector_AS( &isn, &ifg, .10, 4, Cases[i].vector_type );
    if( ifg.clear_rng.bgn < Cases[i].answer.bgn - 4 ||
        ifg.clear_rng.end > Cases[i].answer.end + 4 )
    {
      // GenericMesg gen;

      fprintf( stderr, "Test %d failed: Expected %d,%d & got %d,%d\n", i + 1,
               Cases[i].answer.bgn, Cases[i].answer.end,
               ifg.clear_rng.bgn, ifg.clear_rng.end );
      fprintf( stderr, "Frag: %s\n", ifg.sequence );
      fprintf( stderr, "Vect: %s\n", isn.sequence );
/*      gen.t = MESG_ISN;
      gen.m = &isn;
      WriteProtoMesg_AS( stderr, &gen );
      gen.t = MESG_IFG;
      gen.m = &ifg;
      WriteProtoMesg_AS( stderr, &gen );
*/
    }
    else
    {
      fprintf( stderr, "Test %d passed: Expected %d,%d & got %d,%d\n", i + 1,
               Cases[i].answer.bgn, Cases[i].answer.end,
               ifg.clear_rng.bgn, ifg.clear_rng.end );
    }
  }      
  }

  return 0;
}





