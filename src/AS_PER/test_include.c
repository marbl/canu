
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


#include <sys/types.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"


typedef struct {
  // BLANK  LINES SHOW GROUPING INTO 8-byte BLOCKS

  uint   deleted:1;
  uint   readType:8;
  uint   hasQuality:1;
  uint   numScreenMatches:16; /* number of screen matches */
  uint   hasModifiedClearRegion:1;  // never used as of Oct 2001 - Jason
#if (FRAGSTORE_VERSION >= VERSION_OF_FRAGSTORE_WITH_MODIFIED_CLEARRANGES )
  uint   hasOVLClearRegion:1; 
  uint   hasCNSClearRegion:1; 
  uint   hasCGWClearRegion:1; 
  uint   spare1:2;
#else
  uint   spare1:5;
#endif
  VLSTRING_SIZE_T clearRegionStart;

  VLSTRING_SIZE_T clearRegionEnd;
#if (FRAGSTORE_VERSION >= VERSION_OF_FRAGSTORE_WITH_MODIFIED_CLEARRANGES )
  VLSTRING_SIZE_T ovlRegionStart; 

  VLSTRING_SIZE_T ovlRegionEnd; 
  VLSTRING_SIZE_T cnsRegionStart; 

  VLSTRING_SIZE_T cnsRegionEnd; 
  VLSTRING_SIZE_T cgwRegionStart; 

  VLSTRING_SIZE_T cgwRegionEnd; 
#endif
  CDS_IID_t readIndex;         /* Internal ID of this read */

  CDS_UID_t accID;             /* Accession ID of this read */

  uint64 sequenceOffset;    /* Offset of the sequence/quality data in the seq Store */

  uint64 sourceOffset;      /* Offset of the source in the source, localePos, and screen Matches  */

#if (FRAGSTORE_VERSION > FRAGSTORE_VERSION_PRODUCTION)
   localIndex;        /* Local ID of this read, assigned by appendFragStore */
#else
  uint32 blankPadTo8byteword;
#endif
  time_t entryTime;

}ShortFragRecord;


void dummy_fcn(void){
  printf("dummy\n");
}
