// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Clark Mobarry
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#include "halign.h"

//static char Usage[] = "%s seq1 seq2";

int main(int argc, char *argv[])
{
   char *seq1, *seq2;
   int   len1, len2;
   int   offset1, offset2;
   H_Alignment_t * aln_ptr = NULL;

   seq1 = "ATCGTCCGGATGAAAATGTCTCGGGGGGGGGGGTCGGG";
   seq2 = "ATCGTCTGGATGAAAAAGTCTCAAGGG";
   len1 = strlen(seq1);
   len2 = strlen(seq2);
   
   // This is the setting for the coordinate system to be just the match.
   offset1 = 0;
   offset2 = 0;

   // Sequence coordinates are base-based, starting from 0
   halignStart( seq1+offset1, // This is the first base in the comparison.
	   seq2+offset2,
	   offset1, offset2,
	   len1, len2,
	   &aln_ptr);
   
#if 0
   printUngappedAlign(aln_ptr);
   printUngappedAlignSharpEnds(aln_ptr);
#endif

   printUngappedAlignSharpEndsOnConsole(aln_ptr, seq1, seq2, 0);
   printUngappedAlignSharpEndsOnConsole(aln_ptr, seq1, seq2, 1);
   printUngappedAlignSharpEndsOnConsole(aln_ptr, seq1, seq2, 2);

   if(1){
     int bgn1, bgn2, len1, len2, nmat;
     printf("Try iterator\n");
     while(iterateUngappedAlignSharpEnds(aln_ptr, bgn1, bgn2, len1, len2, nmat)) {
       printf("%d %d %d %d\n", bgn1, bgn2, len1, len2 );
     }
     printf("Tried iterator\n");
   }

   if(aln_ptr != NULL) Free_align(aln_ptr);
   // Must call for each halign() but after printing output.
   
   exit(0);
}
