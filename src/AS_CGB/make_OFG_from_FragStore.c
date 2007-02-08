
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
#include <stdlib.h> 
#include <stdio.h> 
#include <assert.h>
#include <unistd.h> /* man 3 getopt */
#include "AS_global.h" 
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"

#define READ_FRAG_SOURCE_FIELD
#define MAX_SOURCE_LENGTH  1000

int main(int argc, char *argv[]){

  CDS_COORD_t begin = -1, end = -1;

  if(argc  < 2 ){
    fprintf(stderr,"Usage: %s [-b <firstElem>] [-e <lastElem>] <StorePath> \n", argv[0]);
    exit(1);
  }

  {
    int ch;
    while ((ch = getopt(argc, argv, "b:e:l")) != EOF){
      switch(ch) {
        case 'e':
          end = atoi(optarg);
          fprintf(stderr,"* end = "F_COORD"\n", end);
          break;
        case 'b':
          begin = atoi(optarg);
          fprintf(stderr,"* begin = "F_COORD"\n", begin);
          break;
        default:
          fprintf(stderr,"* Unknown option %s\n", optarg);
          break;
      }
    }
  }

  char *frag_store_name = argv[optind];

  FragStoreHandle source;

  source = openFragStore(frag_store_name,"r");

  if(source == NULLSTOREHANDLE){
    exit(1);
  }
  
  if(begin < getFirstElemFragStore(source) || begin > getLastElemFragStore(source)){
    begin = getFirstElemFragStore(source);
  }
  if(end > getLastElemFragStore(source) || end < getFirstElemFragStore(source)){
    end = getLastElemFragStore(source);
  }

  {
    ReadStructp myRead = new_ReadStruct();
    Fragment_ID    uid;
    IntFragment_ID iid;
    FragType read_type;
    uint32   clr_rng_bgn, clr_rng_end, clr_rng_len;  // Clear range
    time_t   etm;
    char     src_text[MAX_SOURCE_LENGTH]= {0};
    int      src_len = MAX_SOURCE_LENGTH;
    
    int iret, ii;
    for(ii = begin; ii <= end; ii++){
      uint32 isDeleted;
      
#ifndef READ_FRAG_SOURCE_FIELD
      getFragStore(source, ii, FRAG_S_FIXED, myRead);
#else
      getFragStore(source, ii, FRAG_S_SOURCE, myRead);
#endif

      iret  = 0;
      iret += getIsDeleted_ReadStruct( myRead, &isDeleted);
      iret += getAccID_ReadStruct( myRead, &uid);
      iret += getReadIndex_ReadStruct( myRead, &iid);
      iret += getReadType_ReadStruct( myRead, &read_type);

#ifdef READ_FRAG_SOURCE_FIELD
      iret += getSource_ReadStruct( myRead, src_text, src_len);
#endif      
      iret += getEntryTime_ReadStruct( myRead, &etm);
      iret += getClearRegion_ReadStruct( myRead, &clr_rng_bgn, &clr_rng_end,
					READSTRUCT_LATEST);  // Clear range

      assert(iret == 0);

      assert(clr_rng_end >= clr_rng_bgn);
      clr_rng_len = clr_rng_end - clr_rng_bgn;


      if (isDeleted) {
        printf("{OFG\nact:D\nacc:("F_UID","F_IID")\n}\n", uid, iid);
      } else { // (!isDeleted)

        switch(read_type) {
          case AS_READ:    //  Celera Read
          case AS_EXTR:    //  External WGS read
          case AS_TRNR:    //  Transposon library read
            printf("{OFG\n"
                   "act:A\n"
                   "acc:("F_UID","F_IID")\n"
                   "typ:%c\n"
                   "src:\n"
                   "%s.\n"
                   "etm:"F_TIME_T"\n"
                   "clr:0,"F_U32"\n"
                   "}\n",
                   uid, iid, read_type, src_text, etm, clr_rng_len );
            break;
          default:
            fprintf(stderr,"Unsupported fragment type ASCII:%c Decimal:%d\n",
                    read_type, read_type);
            assert(FALSE);
            break;
        }
      }
    }
  }

  exit(0);
}
