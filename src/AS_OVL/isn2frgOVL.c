
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
/*************************************************
* Module:  isn2frg.c
* Description:
*   Reads a stream of screen items from standard input
*   and converts them to fragment messages.
*
*    Programmer:  A. Delcher
*       Written:  23 Mar 2000
*  Last Revised:  
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: isn2frgOVL.c,v 1.5 2005-09-15 15:20:16 eliv Exp $
 * $Revision: 1.5 $
*/

static char fileID[] = "$Id: isn2frgOVL.c,v 1.5 2005-09-15 15:20:16 eliv Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include  "AS_OVL_delcher.h"


char  Filter
    (char ch);



int main  (int argc, char * argv [])

  {
   FragMesg  frag_msg;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   MesgWriter  write_msg_fn;
   GenericMesg  * pmesg;
   AuditLine  audit_line;
   AuditMesg  * new_adt_mesg;
   char  label_line [1000];
   int  i;
   char  qual [AS_READ_MAX_LEN + 1];


   read_msg_fn = (MesgReader)InputFileType_AS (stdin);
   write_msg_fn = (MesgWriter)OutputFileType_AS (AS_PROTO_OUTPUT);

   pmesg = (GenericMesg *) Safe_malloc (sizeof (GenericMesg));
   pmesg -> t = MESG_ADT;
   pmesg -> m = (AuditMesg *) Safe_malloc (sizeof (AuditMesg));
   new_adt_mesg = pmesg -> m;
   new_adt_mesg -> list = & audit_line;
      

   while  (read_msg_fn (stdin, & gmesg) != EOF && gmesg != NULL)
     switch  (gmesg -> t)
       {
        case  MESG_ADT :
          {
           AuditMesg  * adt_mesg = gmesg -> m;

           sprintf (label_line, "%s", argv [0]);
           AppendAuditLine_AS (adt_mesg, & audit_line, time (0), "get-subgraph",
                               "$Revision: 1.5 $", label_line);
           write_msg_fn (stdout, gmesg);
           break;
          }

        case  MESG_ISN :
          {
           InternalScreenItemMesg  * isn_mesg = gmesg -> m;
           int  len;
          
           len = strlen (isn_mesg -> sequence);
           assert (len <= AS_READ_MAX_LEN);

           for  (i = 0;  i < len;  i ++)
             {
              isn_mesg -> sequence [i] = Filter (isn_mesg -> sequence [i]);
              qual [i] = 'N';
             }
           qual [i] = '\0';

           frag_msg . action = AS_ADD;
           frag_msg . eaccession = isn_mesg -> eaccession;
           frag_msg . type = AS_READ;
           frag_msg . entry_time = time (NULL);
           frag_msg . clear_rng . bgn = 0;
           frag_msg . clear_rng . end = strlen (isn_mesg -> sequence);
           frag_msg . source = isn_mesg -> source;
           frag_msg . sequence = isn_mesg -> sequence;
           frag_msg . quality = qual;

           gmesg -> t = MESG_FRG;
           gmesg -> m = & frag_msg;
           write_msg_fn (stdout, gmesg);

           break;
          }

        default :
          write_msg_fn (stdout, gmesg);
       }

   return  0;
  }



char  Filter
    (char ch)

//  Return a single  a, c, g or t  for  ch .

  {
   switch  (tolower (ch))
     {
      case  'a' :
      case  'c' :
      case  'g' :
      case  't' :
        return  ch;
      case  'r' :     // a or g
        return  'g';
      case  'y' :     // c or t
        return  'c';
      case  's' :     // c or g
        return  'c';
      case  'w' :     // a or t
        return  't';
      case  'm' :     // a or c
        return  'c';
      case  'k' :     // g or t
        return  't';
      case  'b' :     // c, g or t
        return  'c';
      case  'd' :     // a, g or t
        return  'g';
      case  'h' :     // a, c or t
        return  'c';
      case  'v' :     // a, c or g
        return  'c';
      default :       // anything
        break;
    }

   return  'a';
  }



