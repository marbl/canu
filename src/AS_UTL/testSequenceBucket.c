
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
#include <assert.h>
#include "AS_UTL_SequenceBucket.h"


int main(void){
  float32 probs[]={0.25,0.25,0.25,0.25};
  
#if 0
  SequenceBucketT *sb = CreateSequenceBucket(2);
  
  seq = "aa";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "ac";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "at";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "ag";
  IncrementBucketPrefix(sb, seq);   
  
  
  seq = "ca";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "cc";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "ct";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "cg";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "ta";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "tc";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "tt";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "tg";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "aa";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "aa";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "aa";
  IncrementBucketPrefix(sb, seq);   
  
  seq = "aa";
  IncrementBucketPrefix(sb, seq);   
  
  for(i = 0; i < sb->numBuckets;i++){
    char buffer[32];
    ComputeBucketString(sb,i,buffer);
    fprintf(stderr,"* Bucket %d for string %s has count %d\n",
            i, buffer, sb->buckets[i]);
  }
  ComputeBucketTheoreticalRates(sb, probs);
  ComputeBucketActualRates(sb);
  CheckBucketRates(sb, 3);
#else
  
  SequenceBucketArrayT *sba = CreateSequenceBucketArray(8);
  
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "aaccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acccttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  IncrementSequenceBucketArrayPrefix(sba, "acctttgg"); 
  
  CheckSequenceBucketArray(sba, probs, 5.0, stderr,"");
  
#endif
  return 0;
}
