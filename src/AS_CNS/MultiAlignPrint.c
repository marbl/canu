
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignStore_CNS.h"
#include "Array_CNS.h"

void 
PrintMultiAlignT(FILE *out,
	         MultiAlignT *ma,
	         GateKeeperStore *frag_store, 
	         int show_qv, 
	         int dots,
                 uint32 clrrng_flag)  {

  char frgTypeDisplay;
  FragType frgTypeData;
  int depth;
  int rc,i;
  int window;
  int length;
  char **multia=NULL; 
  int **idarray;
  int **oriarray;
  char *consensus = Getchar(ma->consensus,0);
  char *quality   = Getchar(ma->quality,0);

  fragRecord rsp;

  length = strlen(consensus);
   
  rc = IMP2Array(GetIntMultiPos(ma->f_list,0),
                 GetNumIntMultiPoss(ma->f_list),
                 GetNumchars(ma->consensus),
                 frag_store,
                 &depth,
                 &multia,
                 &idarray,
                 &oriarray,
                 0,
                 clrrng_flag);

  if (rc) {
    fprintf(out,"<<< begin Contig %d >>>",ma->maID);;

    int ungapped=1;
    int tick=1;

    for (window=0;window<length;) {
      int row_id=0;
      int rowind=0;
      int orient=0;
      int rowlen=0;
      char *rowchar=consensus+window;

      fprintf(out, "\n");
      fprintf(out, "\n");
      fprintf(out, "<<<  Contig %d, gapped length: %d  >>>\n",ma->maID, length);
      fprintf(out, "%d gapped\n",window+1);
      fprintf(out, "         |         |         |         |         |         |         |         |         |         |\n");
      fprintf(out, "%d ungapped\n",ungapped+tick-1);

      rowlen = (window+100 < length)?100:length-window;

      for (rowind=0;rowind<rowlen;rowind++,rowchar++){
        if ( tick==10 ) {
          ungapped+=10;
          tick=0;
        }
        if ( tick==0 && *rowchar!='-') {
          fprintf(out,"u");
        } else {
          fprintf(out," ");
        }
        if (*rowchar!='-') {
          tick++;
        } 
      }     

      fprintf(out,"\n");
      fprintf(out,"%-100.100s  cns  (uid,iid) type\n", consensus+window);
      if (show_qv)
        fprintf(out,"%-100.100s  qlt\n", quality+window);

      fprintf(out,"___________________________________________________________________________________________________________________________________\n");

      for (i=0;i<depth;i++) {
        int j;

        if (multia[2*i] == NULL)
          continue;

        char *nonblank = strpbrk(multia[2*i]+window,"ACGT");
        if (nonblank == NULL || nonblank - (multia[2*i]+window) > 100)
          continue;

        for (j=0;j<100;j++) {
          if ( window+j> length) break;
          if ( dots && *(multia[2*i]+window+j) == *(consensus+window+j) ) {
            *(multia[2*i]+window+j) = '.';
            *(multia[2*i+1]+window+j) = ' ';
          } else {
            *(multia[2*i]+window+j) = tolower(*(multia[2*i]+window+j));
          }
        }
        
        {
          int last = (window+99< length-1)?window+99:length-1;
          if ( *(idarray[i]+last) == 0 ) {
            row_id = *(idarray[i]+window);
            orient = *(oriarray[i]+window);
          } else {
            row_id = *(idarray[i]+last);
            orient = *(oriarray[i]+last);
          }
        }

        // Look up UID for row_id
        if ( row_id > 0 ) {
          getFrag(frag_store, row_id, &rsp, FRAG_S_INF);

          //getReadType_ReadStruct(rsp, &frgTypeData);
          //
          // AS_READ is normally 'R' but for reports, we were asked
          // to show ' ' for these perponderant normal reads.
          //
          //  placeholder for future expansion to different types...
          //
          frgTypeData = AS_READ;
          frgTypeDisplay = ' ';

          fprintf(out, "%-100.100s   %c   (%s,%d) %c\n",
                  multia[2*i]+window,
                  (orient>0)?'>':'<',
                  AS_UID_toString(getFragRecordUID(&rsp)),
                  row_id,
                  frgTypeDisplay);

          if (show_qv)
            fprintf(out, "%-100.100s   %c   (%s,%d) %c\n",
                    multia[2*i+1]+window,
                    (orient>0)?'>':'<',
                    AS_UID_toString(getFragRecordUID(&rsp)),
                    row_id,
                    frgTypeDisplay);
        }
      }
      window+=100;
    }
    fprintf(out,"\n<<< end Contig %d >>>\n", ma->maID);
  } else {
    fprintf(stderr,"Error returned from MultiAlignT2Array.\n");
  }

  if (multia) {
    for (i=0;i<2*depth;i++)
      safe_free(multia[i]);

    safe_free(multia);

    for (i=0;i<depth;i++) {
      safe_free(idarray[i]);
      safe_free(oriarray[i]);
    }

    safe_free(idarray);
    safe_free(oriarray);
  }
}
