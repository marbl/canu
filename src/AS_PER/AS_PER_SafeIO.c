
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
/*************************************************************************
 Module:  AS_PER_SafeIO
 Description:
   Safe Read and Write operations.

 Assumptions:
      

 *************************************************************************/

/* RCS Info
 * $Date: 2004-04-14 13:52:55 $
 * $Id: AS_PER_SafeIO.c,v 1.1.1.1 2004-04-14 13:52:55 catmandew Exp $
 * $Revision: 1.1.1.1 $
 *
 */

#include <stdlib.h>
#include <stdio.h>

#include "AS_global.h"

#include "AS_PER_SafeIO.h"
//#define DEBUG 

int safeRead(FILE *fp, void *p, size_t readSize){
  size_t bytesRead = 0;
  size_t readResult;

  while(bytesRead < readSize){

    clearerr(fp);
    readResult = fread((void *)((char *)p+ bytesRead),1,readSize-bytesRead,fp);
    if(feof(fp) || readResult == 0){
#if 1
      fprintf(stderr," *** EOF encountered in safeRead read " F_SIZE_T "/" F_SIZE_T "...at offset " F_OFF_T " returning TRUE\n",
	      bytesRead, readSize, CDS_FTELL(fp));
      //      exit(1);
#endif
      return TRUE;
    }else if(ferror(fp)){
      fprintf(stderr," *** Error encountered in safeRead...exiting\n");
      exit(1);
    }
#ifdef DEBUG
        fprintf(stderr," safeRead(" F_SIZE_T ") read " F_SIZE_T "/" F_SIZE_T " bytes offset = " F_OFF_T "\n",
	  readResult, bytesRead, readSize, CDS_FTELL(fp));
#endif
    bytesRead += readResult;
  }
  return(0);

}


int safeWrite(FILE *fp, void *p, size_t writeSize){
  size_t bytesWritten = 0;
  size_t writeResult;
#ifdef DEBUG
  off_t previousOffset;
#endif

  while(bytesWritten < writeSize){

    clearerr(fp);
#ifdef DEBUG
    previousOffset = CDS_FTELL(fp);
#endif
    writeResult = fwrite((void *)((char *)p+ bytesWritten),1,
                         writeSize-bytesWritten,fp);
#ifdef DEBUG
    fprintf(stderr," safeWrite(" F_SIZE_T ") wrote " F_SIZE_T " bytes at offset %ld, now at " F_OFF_T "\n",
	    writeResult, writeSize, previousOffset, CDS_FTELL(fp));
#endif
    if(ferror(fp)){
      char buffer[256];
      sprintf(buffer,"***safeWrite: fileno = %d error = %d: ",
	      fileno(fp), ferror(fp));
      perror(buffer);
      
      exit(1);
    }
    bytesWritten += writeResult;
  }
  return(0);
}
