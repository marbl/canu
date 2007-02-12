/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2005, J. Craig Venter Institute
 * Author: Brian Walenz
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
#include <unistd.h>

#include "AS_global.h" 
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"

#define MAX_SOURCE_LENGTH  1000

void
reapFragStore(int32 begin, int32 end, char *gkpname) {  
  GateKeeperStore *source;
  ReadStruct      *myRead;
  int              ii;

  source = openGateKeeperStore(gkpname, FALSE);
  if (source == NULL) {
    fprintf(stderr, "Failed to open gatekeeper store (\"%s\")\n", gkpname);
    exit(1);
  }

  myRead = new_ReadStruct();

  if ((begin < getFirstElemFragStore(source)) ||
      (begin > getLastElemFragStore(source)))
    begin = getFirstElemFragStore(source);
  if ((end > getLastElemFragStore(source)) ||
      (end < getFirstElemFragStore(source)))
    end = getLastElemFragStore(source);

  for (ii = begin; ii <= end; ii++) {
    uint32           isDeleted = 0;
    Fragment_ID      uid = 0;
    IntFragment_ID   iid = 0;
    FragType         read_type = AS_READ;
    time_t           etm = 0;
    uint32           clr_rng_bgn = 0;
    uint32           clr_rng_end = 0;
    int              iret = 0;
      
    getFrag(source, ii, myRead, FRAG_S_INF);

    iret += getIsDeleted_ReadStruct(myRead, &isDeleted);
    iret += getAccID_ReadStruct(myRead, &uid);
    iret += getReadIndex_ReadStruct(myRead, &iid);
    iret += getClearRegion_ReadStruct(myRead, &clr_rng_bgn, &clr_rng_end, READSTRUCT_LATEST);

    //  We don't expect any of these to fail, so lump them all together.
    if (iret)
      fprintf(stderr, "Failed to *_ReadStruct %d times.\n", iret), exit(1);

    //  This one is easier to fail, especially if we are playing with the clear ranges.
    //
    if (clr_rng_bgn > clr_rng_end)
      fprintf(stderr, "Clear range is invalid for ("F_UID","F_IID"): begin="F_U32" > end="F_U32"\n", uid, iid, clr_rng_bgn, clr_rng_end), exit(1);

    if (isDeleted) {
      printf ("{OFG\n"
              "act:D\n"
              "acc:("F_UID","F_IID")\n"
              "}\n",
              uid, iid);
    } else {
      //printf(" "F_IID" %c "F_U32" "F_U32" "F_U32" "F_U32"\n",
      //       iid, read_type, isDeleted, clr_rng_bgn, clr_rng_end, clr_rng_len);

      // definitions in cds/AS/src/AS_MSG/AS_MSG_pmesg.h

      switch (read_type) {
        case AS_READ:   // Celera Read
        case AS_EXTR:   // External WGS read
        case AS_TRNR:   // Transposon library read
          printf("{OFG\n"
                 "act:A\n"
                 "acc:("F_UID","F_IID")\n"
                 "typ:%c\n"
                 "src:\n"
                 ".\n"
                 "etm:"F_TIME_T"\n"
                 "clr:"F_U32","F_U32"\n"
                 "}\n",
                 uid, iid, read_type, etm, clr_rng_bgn, clr_rng_end);
          break;
        default:
          fprintf(stderr,"Unsupported fragment type ASCII:%c Decimal:%d\n", read_type, read_type);
          exit(1);
          break;
      }
    }
  }
}


int main(int argc, char *argv[]) {
  CDS_COORD_t begin = -1;
  CDS_COORD_t end   = -1;
  int ch;

  if (argc < 2) {
    fprintf(stderr,"Usage: %s [-b <firstElem>] [-e <lastElem>] <StorePath> \n", argv[0]);
    exit(1);
  }

  while ((ch = getopt(argc, argv, "b:e:l")) != EOF) {
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

  reapFragStore(begin, end, argv[optind]);

  exit(0);
}





