
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
/* ScaffoldMap
 *    Prints an iid-uid map of the SCF messages in a protoIO file
 *
 * $Id: ScaffoldMap.c,v 1.6 2006-09-21 21:34:00 brianwalenz Exp $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "AS_global.h"
#include <assert.h>

#define MAX_MESG (NUM_OF_REC_TYPES)
int main(int argc, char *argv[])
{ 
  GenericMesg *pmesg;
  int i = 0;
  MesgReader reader = (MesgReader)InputFileType_AS(stdin);

  while (reader(stdin,&pmesg) != EOF){
    assert(pmesg->t <= MAX_MESG );
    switch(pmesg->t){
      default:
        break;
      case MESG_SCF:
        {
          SnapScaffoldMesg *m = pmesg->m;
          if(i++ == 0){
            fprintf(stdout,"* Int  Ext IDs of Scaffolds\n");
          }
          fprintf(stdout,F_IID " " F_UID "\n",
                  m->iaccession, m->eaccession);
        }
    }
  }
  exit (0);
}
