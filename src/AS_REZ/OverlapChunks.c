
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
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"

//
// AS_CGW
//
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ChiSquareTest_CGW.h"

//
// AS_REZ
//
#include "DataTypesREZ.h"
#include "UtilsREZ.h"
#include "CommonREZ.h"
#include "SubgraphREZ.h"
#include "GapWalkerREZ.h"
#include "FbacREZ.h"
#include "dpc_CNS.h"
//
// vars
//

int TestOverlapChunks(void);

Overlap* NewerOverlapChunks( char *seq1, char *seq2,
							 ChunkOrientationType orientation, 
							 int beg, int end)
{
  Overlap *omesg;
  double erate, thresh;
  int minlen;
  int flip = 0;
  // int length1;
  // char *seq1_orig;
  
  erate = CGW_DP_ERATE;
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN;
  
  // if the orientation is BA_AB or BA_BA, we need to reverse complement the first contig
  if (orientation == BA_AB || orientation == BA_BA)
  {
	// length1 = strlen( seq1 );
	// seq1_orig = (char *) malloc( length1 + 1 );
	// strcpy (seq1_orig, seq1);
	Complement_Seq( seq1 );
  }
  
  // if the orientation is AB_BA or BA_BA, we need to set the flip variable for the second contig
  if (orientation == AB_BA || orientation == BA_BA)
	flip = 1;

  // beg and end are essentially bounds on the a-hang
  omesg = DP_Compare(seq1, seq2,
					 beg, end, (int) flip,
					 erate, thresh, minlen,
					 AS_FIND_ALIGN);

  // be nice to whoever called this routine, return seq1 to its original state
  if (orientation == BA_AB || orientation == BA_BA)
	Complement_Seq( seq1 );
  
  // omesg->begpos is the a-hang, omesg->endpos is the b-hang
  return omesg;
}

int TestOverlapChunks()
{
  Overlap *omesg;
  double erate, thresh;
  int minlen;
  int min_ahang, max_ahang;
  ChunkOrientationType orientation;
  char *test1, *test2;
  int len1, len2;
#if 0
  int flip = 0;
#endif

#if 0
  // AB_AB begpos -9 endpos 22
  char *consensus1 =          "aaaaacccccccccccccccccggggggggggggggggggggggggttttttttttttttttttttttttt";
  char *consensus2 = "aaaaaaaaaaaaaacccccccccccccccccggggggggggggggggggggggggtttttttttttttttttttttttttaaaaaaaaaaaaaaaaaaaaaa";
  orientation = AB_AB;
  min_ahang = -7;
  max_ahang = -5;
#endif

#if 1
  // AB_AB begpos -9 endpos 22
  char *consensus1 =          "aaaaacccccccccccccccccggggggggggggggggggggggggttttttttttttttttttttttttt";
  char *consensus2 = "aaaaaaaaaaaaaacccccccccccccccccggggggggggggggggggggggggtttttttttttttttttttttttttaaaaaaaaaaaaaaaaaaaaaa";
  orientation = AB_AB;
  min_ahang = -7;
  max_ahang = -5;
#endif

#if 0
  // AB_BA begpos -9 endpos 22
  char *consensus1 =                          "aaaaacccccccccccccccccggggggggggggggggggggggggttttttttttttttttttttttttt";
  char *consensus2 = "ttttttttttttttttttttttaaaaaaaaaaaaaaaaaaaaaaaaaccccccccccccccccccccccccgggggggggggggggggtttttttttttttt";
  orientation = AB_BA;
  min_ahang = -7;
  max_ahang = -5;
  flip = 1;
#endif

#if 0
  // BA_AB
  // char *consensus1 =       "aaaaaaaaaaaaaaaaaaaaaaaaaccccccccccccccccccccccccgggggggggggggggggttttt";
  // char *consensus2 = "aaaaaaaaaaaaaacccccccccccccccccggggggggggggggggggggggggtttttttttttttttttttttttttaaaaaaaaaaaaaaaaaaaaaa";
  orientation = BA_AB;
#endif

#if 0
  // BA_BA
  // char *consensus1 =                       "aaaaaaaaaaaaaaaaaaaaaaaaaccccccccccccccccccccccccgggggggggggggggggttttt";
  // char *consensus2 = "ttttttttttttttttttttttaaaaaaaaaaaaaaaaaaaaaaaaaccccccccccccccccccccccccgggggggggggggggggtttttttttttttt";
  char *consensus1 =          "aaaaacccccccccccccccccggggggggggggggggggggggggttttttttttttttttttttttttt";
  char *consensus2 = "aaaaaaaaaaaaaacccccccccccccccccggggggggggggggggggggggggtttttttttttttttttttttttttaaaaaaaaaaaaaaaaaaaaaa";
  orientation = BA_BA;
  min_ahang = -7;
  max_ahang = -5;
  flip = 0;
#endif

  len1 = strlen(consensus1);
  test1 = (char *) malloc(len1);
  strncpy( test1, consensus1, len1);
  
  len2 = strlen(consensus2);
  test2 = (char *) malloc(len2);
  strncpy( test2, consensus2, len2);
  
  erate = CGW_DP_ERATE;
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN;
  
  // omesg = DP_Compare(consensus1, consensus2,
  //omesg = DP_Compare(test1, test2,
  //				 min_ahang, max_ahang, (int) flip,
  //				 erate, thresh, minlen,
  //				 AS_FIND_ALIGN);
  
  omesg = NewerOverlapChunks( test1, test2,
							  orientation,
							  min_ahang, max_ahang);

  if (omesg != NULL)
  {
	fprintf( stderr, "omesg->begpos = %d, omesg->endpos = %d\n", omesg->begpos, omesg->endpos);
	fprintf( stderr, "omesg->length = %d\n", omesg->length);
	fprintf( stderr, "omesg->diffs = %d\n", omesg->diffs);
	fprintf( stderr, "omesg->comp = %d\n", omesg->comp);
	return 0;
  }
  else
	fprintf( stderr, "No alignment found!\n");
  return 1;
}

int main(void)
{
  TestOverlapChunks();
  return 0;
}

