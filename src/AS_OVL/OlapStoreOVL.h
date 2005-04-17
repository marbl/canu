
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
/* RCS info
 * $Id: OlapStoreOVL.h,v 1.5 2005-04-17 20:11:58 ahalpern Exp $
*/


#define  OVL_STORE_VERSION       0

#ifndef  __OLAPSTOREOVL_H_INCLUDED
#define  __OLAPSTOREOVL_H_INCLUDED


#include <sys/param.h>
//#include  "AS_OVL_delcher.h"
#include  "AutoScreenOVL.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_fragStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_UTL_version.h"

#define  AIID_BITS               31
#define  AHANG_BITS              16
#define  ERATE_BITS              16
    //  Number of bits to store integer versions of error rates
#define  MAX_ERATE               ((1 << ERATE_BITS) - 1)
    //  Maximum value allowed for integer versions of error rates
#define  MAX_FILENAME_LEN      2000
  // Most characters allowed in a filename
#define  OFFSET_FILENAME       "/offset.olap"
  // Name of file containing offsets in an overlap store

#define AIID_TYPE unsigned
#define AHANG_TYPE signed int
#define ERATE_TYPE unsigned

#ifndef BYTE_ORDER
# error BYTE_ORDER macro was not defined by any include file
#endif

typedef  union
  {
   uint32  header;
   struct
     {
#if BYTE_ORDER == LITTLE_ENDIAN
      unsigned int  record_size : 16;
      unsigned int  version : 16;
#else
      unsigned int  version : 16;
      unsigned int  record_size : 16;
#endif
     }  tag;
  }  OVL_Store_ID_t;

typedef  struct
  {

   AIID_TYPE  a_iid : AIID_BITS;
   unsigned  unused : 1;
   AIID_TYPE  b_iid : AIID_BITS;
   unsigned  flipped : 1;
   AHANG_TYPE  a_hang : AHANG_BITS;
   AHANG_TYPE  b_hang : AHANG_BITS;
#if  OVL_STORE_VERSION == 1
   AHANG_TYPE  corr_a_hang : AHANG_BITS;
   AHANG_TYPE  corr_b_hang : AHANG_BITS;
#endif
   ERATE_TYPE  orig_erate : ERATE_BITS;   // original error rate (aka quality)
   ERATE_TYPE  corr_erate : ERATE_BITS;   // error rate after fragment correction
  }  Long_Olap_Data_t;

typedef  struct
  {
   AIID_TYPE  b_iid : AIID_BITS;
   unsigned  flipped : 1;
   AHANG_TYPE  a_hang : AHANG_BITS;
   AHANG_TYPE  b_hang : AHANG_BITS;
#if  OVL_STORE_VERSION == 1
   AHANG_TYPE  corr_a_hang : AHANG_BITS;
   AHANG_TYPE  corr_b_hang : AHANG_BITS;
#endif
   ERATE_TYPE  orig_erate : ERATE_BITS;   // original error rate (aka quality)
   ERATE_TYPE  corr_erate : ERATE_BITS;   // error rate after fragment correction
  }  Short_Olap_Data_t;

typedef  struct
  {
   char  * name;
   FILE  * offset_fp;
   uint32  max_frag;
   uint32  frags_per_file;
  }  OVL_Store_t;

typedef  struct
  {
   OVL_Store_t  * store;
   FILE  * fp;
   uint32  start_id, stop_id, curr_id;
   int  curr_file_index;
   uint32  curr_offset;
   uint32  * offset_buff;
  }  OVL_Stream_t;


float32  Expand_Quality
    (int q);
void  Free_OVL_Store
    (OVL_Store_t * store);
void  Free_OVL_Stream
    (OVL_Stream_t * stream);
int16  Get_Int_Quality
    (int q);
void  Init_OVL_Stream
    (OVL_Stream_t * stream, uint32 first, uint32 last, OVL_Store_t * store);
void  Init_OVL_Stream_Intra_Frg
    (OVL_Stream_t * stream, uint32 iid, int skipped_ovls, OVL_Store_t * store);
uint32  Last_Frag_In_OVL_Store
    (OVL_Store_t * store);
uint32 *  Load_Frag_Offsets
    (OVL_Store_t * store);
FILE *  Local_File_Open
    (const char * filename, const char * mode);
OVL_Store_t *  New_OVL_Store
    (void);
OVL_Stream_t *  New_OVL_Stream
    (void);
void  Renew_OVL_Stream
    (OVL_Stream_t *stream);
int  Next_From_OVL_Stream
    (Long_Olap_Data_t * olap, OVL_Stream_t * stream);
int  Open_OVL_Store
    (OVL_Store_t * store, const char * path);
int  Shrink_Quality
    (float32 q);

#endif
