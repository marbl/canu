
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

/**********************************************************************
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDclient_stand_alone.c,v $
$Revision: 1.2 $
$Date: 2005-03-22 19:08:55 $
$Name: not supported by cvs2svn $
$Author: jason_miller $
$Log: not supported by cvs2svn $
Revision 1.3  2004/09/10 12:31:43  mschatz
Add standard copyright notice

Revision 1.2  2004/09/09 22:38:28  mschatz
Type safety

Revision 1.1  2004/06/24 12:51:06  mpop
Added AS_UID

Revision 1.2  2003/05/09 21:04:03  mpop
Dos2unixed all files.
Modified c_make.as to set SEP_PATH relative to LOCAL_WORK

Revision 1.1.1.1  2003/05/08 18:40:12  aaronhalpern
versions from TIGR

Revision 1.2  2001/09/25 23:03:20  mpop
Dos2Unixed

Revision 1.1.1.1  2001/09/25 20:21:05  mpop
Celera Assembler

Revision 1.1  1999/10/13 17:19:10  sdmurphy
baseline


**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "cds.h"
#include "SYS_UIDcommon.h"
#include "SYS_UIDclient.h"

cds_int32 main(cds_int32 argc, char** argv)
{
   cds_uint64             uid_interval[4];
   cds_int32              uid_status;
   cds_int32              msec_delay = 0;
   cds_int32              count = 0;
   cds_int32              print_flag = 0;
   cds_int32              i;
   cds_uint64             block_size;
   cds_uint64             a_new_uid;
   cds_int32              a_new_uid_status;
   cds_int32              j;
   char               err_str[300];
   cds_uint64             max_block_size = 0L;

   // check
   if (argc < 2)
   {
      sprintf(err_str, "usage: %s <block_size> \n", argv[0]);
      fprintf(stderr, err_str);
      exit(0);
   }
   
   // parse
   block_size      = strtoul(argv[1], (char**)NULL, 10);

   uid_status = SYS_UIDgetMaxUIDSize(&max_block_size);
   if (max_block_size < block_size) {
     uid_status = UID_CODE_BLOCK_TOO_LARGE;
   }
   if (uid_status != UID_CODE_OK) {
     printf(F_S32" 0 0 0 0\n", uid_status);
     exit(1);
   }

   SYS_UIDsetUIDSize(block_size);
   uid_status = SYS_UIDgetNewUIDInterval(uid_interval);
   printf(F_S32" "F_S64" "F_S64" "F_S64" "F_S64"\n",
		uid_status, 
		uid_interval[0], 
		uid_interval[1],
		uid_interval[2], 
		uid_interval[3]);
   if (uid_status == UID_CODE_OK) {
     exit(0);
   } else {
     exit(1);
   }
}


