
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
#include <stdlib.h>

#include "SYS_UIDcommon.h"
#include "SYS_UIDclient.h"

cds_int32 main(cds_int32 argc, char** argv)
{
   cds_uint64             uid_interval[4];
   cds_int32              uid_status;
   cds_int32              msec_delay = 0;
   cds_int32              count = 0;
   cds_int32              print_flag = 0;
   cds_uint32             us;
   cds_int32              i;
   cds_uint64             block_size;
   cds_uint64             a_new_uid;
   cds_int32              a_new_uid_status;
   cds_int32              j;
   char               err_str[300];
   cds_uint64             max_block_size = 0L;

   // check
   if (argc < 3)
   {
      sprintf(err_str,
	      "usage: %s <block_size> <count> <ms_delay> [-p to print]\n",
	      argv[0]);
      fprintf(stderr, err_str);
      return 0;
   }
   
   // parse
   block_size      = strtoul(argv[1], (char**)NULL, 10);
   count           = atoi(argv[2]);
   msec_delay      = atoi(argv[3]);

   if (argc >= 5)
      if (strcmp("-p", argv[4]) == 0)
         print_flag = 1;

   uid_status = SYS_UIDgetMaxUIDSize(&max_block_size);

   if (print_flag)
     printf("max UID status : "F_S32" max size: "F_S64" \n", uid_status, max_block_size);

   // begin loop
   SYS_UIDsetUIDSize(block_size);
   for (i=0;i<count;i++)
   {
      uid_status = SYS_UIDgetNewUIDInterval(uid_interval);
      if (print_flag)
         printf(F_S32": " F_S64" "F_S64" "F_S64" "F_S64"\n", 
		uid_status, 
		uid_interval[0], 
		uid_interval[1],
		uid_interval[2], 
		uid_interval[3]);  

      for (j=0; j<block_size; j++)
      {
         a_new_uid_status = SYS_UIDgetNextUID(&a_new_uid);
         if (print_flag)
	   printf(F_S32"_"F_S64" ", a_new_uid_status, a_new_uid); 
      }
      if (print_flag)
	printf("\n");
      us = msec_delay * 1000;
      usleep(us);
   }
   return 0;
}


