
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

void reapFragStore(int load, int32 begin, int32 end, char * frag_store_name)
{  
  FragStoreHandle source;

  if(load){
    fprintf(stderr,"* LOADing FragStore %s\n", frag_store_name);
    source = loadFragStorePartial(frag_store_name,STREAM_FROMSTART, STREAM_UNTILEND);
  }else{
    fprintf(stderr,"* Opening FragStore %s\n", frag_store_name);
    source = openFragStore(frag_store_name,"r");
  }

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
    Locale_ID     locID;
    IntLocale_ID  lid;
    uint32 locale_pos_bgn, locale_pos_end;
    
    
    int iret, ii;
    for(ii = begin; ii <= end; ii++){
      uint32 isDeleted;
      
#ifndef READ_FRAG_SOURCE_FIELD
      getFragStore(source, ii, FRAG_S_FIXED, myRead);
#else // READ_FRAG_SOURCE_FIELD
      getFragStore(source, ii, FRAG_S_SOURCE, myRead);
#endif // READ_FRAG_SOURCE_FIELD
      // Options: FRAG_S_FIXED, FRAG_S_SOURCE, FRAG_S_SEQUENCE, FRAG_S_ALL
      // dump_ReadStruct(myRead, stdout, FALSE);

      iret = getIsDeleted_ReadStruct( myRead, &isDeleted);
      assert(iret == 0);

      iret = getAccID_ReadStruct( myRead, &uid);
      assert(iret == 0);
      iret = getReadIndex_ReadStruct( myRead, &iid);
      assert(iret == 0);
      iret = getReadType_ReadStruct( myRead, &read_type);
      assert(iret == 0);
      iret = getLocID_ReadStruct( myRead, &locID);
      assert(iret == 0);
#if 0
      iret = getLocalIndex_ReadStruct( myRead, &lid);
#else
      lid = 1;
#endif      
      assert(iret == 0);
      iret = getLocalePos_ReadStruct( myRead, &locale_pos_bgn, &locale_pos_end);
      assert(iret == 0);
#ifdef READ_FRAG_SOURCE_FIELD
      iret = getSource_ReadStruct( myRead, src_text, src_len);
      assert(iret == 0);
#endif      
      iret = getEntryTime_ReadStruct( myRead, &etm);
      assert(iret == 0);
      iret = getClearRegion_ReadStruct( myRead, &clr_rng_bgn, &clr_rng_end,
					READSTRUCT_LATEST);  // Clear range
      assert(iret == 0);
      assert(clr_rng_end >= clr_rng_bgn);
      clr_rng_len = clr_rng_end - clr_rng_bgn;


      if( isDeleted ) {
        printf("{OFG\nact:D\nacc:(" F_UID "," F_IID ")\n}\n",
               uid, iid);
      } else { // (!isDeleted)
#if 0
        printf(" " F_IID " %c " F_U32 " " F_U32 " " F_U32 " " F_U32 "\n",
               iid, read_type, isDeleted, clr_rng_bgn, clr_rng_end, clr_rng_len);
#else
      switch(read_type) {
        // definitions in cds/AS/src/AS_MSG/AS_MSG_pmesg.h

      case AS_READ:  // (int)'R',  // Celera Read
      case AS_EXTR:  // (int)'X',  // External WGS read
      case AS_TRNR:  // (int)'T',  // Transposon library read
      case AS_B_READ: // (int)'G',  // BGLII read
        printf("{OFG\nact:A\nacc:(" F_UID "," F_IID ")\ntyp:%c\nsrc:\n%s.\netm:" F_TIME_T "\nclr:0," F_U32 "\nscn:\n.\n}\n",
               uid, iid, read_type, src_text, etm, clr_rng_len );
        break;
      case AS_EBAC: // (int)'E',  //End of BAC
        printf("{OFG\nact:A\nacc:(" F_UID "," F_IID ")\ntyp:%c\n"
               "loc:(" F_UID "," F_IID ")\n"
               "src:\n%s.\netm:" F_TIME_T "\nclr:0," F_U32 "\nscn:\n.\n}\n",
               uid, iid, read_type,
               locID, lid,
               src_text, etm, clr_rng_len );
        break;
      case AS_FBAC: // (int)'F',  //Finished
        printf("{OFG\nact:A\nacc:(" F_UID "," F_IID ")\ntyp:%c\n"
               "loc:(" F_UID "," F_IID ")\n"
               "sid:(%d,%d)\n"
               "pos:%d,%d\n"
               "src:\n%s.\netm:" F_TIME_T "\nclr:0," F_U32 "\nscn:\n.\n}\n",
               uid, iid, read_type,
               locID, lid,
               0,0, 0,0, // These are incorrect values.
               src_text, etm, clr_rng_len );
        break;
      case AS_UBAC: // (int)'U',  //Unfinished BACs
        printf("{OFG\nact:A\nacc:(" F_UID "," F_IID ")\ntyp:%c\n"
               "loc:(" F_UID "," F_IID ")\n"
               "sid:(%d,%d)\n"
               "btd:(%d,%d)\n"
               "pos:%d,%d\n"
               "src:\n%s.\netm:" F_TIME_T "\nclr:0," F_U32 "\nscn:\n.\n}\n",
               uid, iid, read_type,
               locID, lid,
               0,0,  // These are incorrect values.
               0,0, 0,0,  // These are incorrect values.
               src_text, etm, clr_rng_len );
        break;
      case AS_LBAC: // (int)'L',  //Lightly shotgunned BACs
      case AS_STS:  // (int)'S',  //Sts
        /* The following are never intended to be for Unitigger input. */
      case AS_UNITIG: // (int)'u',  //Unitig
      case AS_CONTIG: // (int)'c',   //Contig
      case AS_BACTIG: // (int) 'B',   // BacTig
      case AS_FULLBAC: // (int)'C'   // Full Bac C = Complete)
      default:
        fprintf(stderr,"Unsupported fragment type ASCII:%c Decimal:%d\n", read_type, read_type);
        assert(FALSE);
        break;
      }
#endif
      }
    }
  }
}


int main(int argc, char *argv[]){

  int load = FALSE;
  CDS_COORD_t begin = -1, end = -1;

  if(argc  < 2 ){
    fprintf(stderr,"Usage: %s [-b <firstElem>] [-e <lastElem>] [-l] <StorePath1> \n",
            argv[0]);
    fprintf(stderr,"   -l option causes frag store to be loaded into memory, rather than opened\n");
    exit(1);
  }

  {
    int ch;
    while ((ch = getopt(argc, argv, "b:e:l")) != EOF){
      switch(ch) {
      case 'l':
        load = TRUE;
        break;
      case 'e':
        end = atoi(optarg);
        fprintf(stderr,"* end = " F_COORD "\n", end);
        break;
      case 'b':
        begin = atoi(optarg);
        fprintf(stderr,"* begin = " F_COORD "\n", begin);
        break;
      default:
        fprintf(stderr,"* Unknown option %s\n", optarg);
        break;
      }
    }
  }

  reapFragStore( load, begin, end, argv[optind]);

  fprintf(stderr,"* Bye Bye\n");
  exit(0);
}





