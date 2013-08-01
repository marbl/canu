
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

static const char *rcsid = "$Id$";

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.H"
#include "AS_UTL_fileIO.H"
#include "AS_UTL_reverseComplement.H"
#include "MultiAlignment_CNS.H"

//  Width of the multialignment display.  100 is convenient for screen display; larger
//  values work well for display in a web browser.
//
uint32 MULTIALIGN_PRINT_WIDTH = 100;

//  Space between fragments on a single line.  If set to a massive value then all
//  fragments appear on different lines.
//
uint32 MULTIALIGN_PRINT_SPACING = 3;



class LaneNode {
public:
  LaneNode() {
    read      = NULL;
    readLen   = 0;
    sequence  = NULL;
    quality   = NULL;
    next      = NULL;
  };

  ~LaneNode() {
    delete [] sequence;
    delete [] quality;
  };

  IntMultiPos      *read;
  int32             readLen;
  char             *sequence;
  char             *quality;

  LaneNode         *next;
};



class Lane {
public:
  Lane() {
    first   = NULL;
    last    = NULL;
    lastcol = 0;
  };

  ~Lane() {
    LaneNode *node = first;
    LaneNode *next = NULL;

    while (node) {
      next = node->next;
      delete node;
      node = next;
    }
  };

  //  Add a node to this lane, if possible.
  bool addNode(LaneNode *node) {
    int32 leftpos = node->read->position.end;
    if (node->read->position.bgn < node->read->position.end)
      leftpos = node->read->position.bgn;

    if ((lastcol > 0) &&
        (leftpos < lastcol + MULTIALIGN_PRINT_SPACING))
      return(false);

    assert(node->next == NULL);

    if (first == NULL) {
      first      = node;
      last       = node;
    } else {
      last->next = node;
      last       = node;
    }

    lastcol = leftpos + node->read->delta_length + node->readLen;

    return(true);
  };

  LaneNode         *first;
  LaneNode         *last;
  int32             lastcol;
};



static
int
IntMultiPositionCmp(const void *l, const void *m) {
  const IntMultiPos *L = (const IntMultiPos *)l;
  const IntMultiPos *M = (const IntMultiPos *)m;

  int32 ltmp = (L->position.bgn < L->position.end) ? L->position.bgn : L->position.end;
  int32 mtmp = (M->position.bgn < M->position.end) ? M->position.bgn : M->position.end;

  if (ltmp == mtmp)
    return 0;

  return((ltmp > mtmp ) ? 1 : -1);
}


#if 0
class alignArray {
public:
  int32    depth;
  char   **seq;
  char   **qlt;
  int32  **IDs;
  int32  **ori;
}
#endif



static
void
IMP2Array(IntMultiPos *imp,
          int32        impLen,
          int32        cnsLen,
          gkStore     *gkp,
          int32       *depth,
          char      ***array,
          int32     ***id_array,
          int32     ***ori_array,
          uint32       clrrng_flag) {

  Lane  *lane     = NULL;
  Lane  *lanes    = new Lane [impLen];
  int32  lanesLen = 0;
  int32  lanesPos = 0;

  // Sort the fragments by leftmost position within contig
  qsort(imp, impLen, sizeof(IntMultiPos), IntMultiPositionCmp);

  //  Load into lanes.

  for (int32 i=0; i<impLen; i++) {
    gkFragment  fsread;
    uint32      clr_bgn;
    uint32      clr_end;
    LaneNode   *node = new LaneNode();

    gkp->gkStore_getFragment(imp[i].ident, &fsread, GKFRAGMENT_QLT);

    fsread.gkFragment_getClearRegion(clr_bgn, clr_end, clrrng_flag);

    //  Hack fix for dumping multialigns when OBT is not run (see tigStore.C:606)
    if ((clr_end < clr_bgn) &&
        (clrrng_flag == AS_READ_CLEAR_OBTCHIMERA))
      fsread.gkFragment_getClearRegion(clr_bgn, clr_end, AS_READ_CLEAR_CLR);

    if (clr_end < clr_bgn)
      fprintf(stderr, "ERROR:  Undefined clear range for fragment %d\n", imp[i].ident), exit(1);

    node->read        = &imp[i];
    node->readLen     = clr_end - clr_bgn;
    node->sequence    = new char [node->readLen + 1];
    node->quality     = new char [node->readLen + 1];

    memcpy(node->sequence, fsread.gkFragment_getSequence() + clr_bgn, sizeof(char) * node->readLen);
    memcpy(node->quality,  fsread.gkFragment_getQuality()  + clr_bgn, sizeof(char) * node->readLen);

    node->sequence[node->readLen] = 0;
    node->quality [node->readLen] = 0;

    if (node->read->position.bgn > node->read->position.end)
      reverseComplement(node->sequence, node->quality, node->readLen);

    //  Try to add this new node to the lanes.  The last iteration will always succeed, adding the
    //  node to a fresh empty lane.

    for (lanesPos=0; lanesPos <= lanesLen; lanesPos++)
      if (lanes[lanesPos].addNode(node))
        break;

    assert(lanesPos <= lanesLen);

    //  If it is added to the last one, increment our cnsLen.

    if (lanesPos == lanesLen)
      lanesLen++;
  }

  //  Process.

  char   **multia = (char **)safe_malloc(2*lanesLen*sizeof(char *));

  for (int32 i=0; i<2*lanesLen; i++) {
    multia[i] = (char *)safe_malloc((cnsLen + 1 + MULTIALIGN_PRINT_WIDTH) * sizeof(char));
    memset(multia[i], ' ', cnsLen + MULTIALIGN_PRINT_WIDTH);
    multia[i][cnsLen] = 0;
  }

  int32  **ia = (int32 **)safe_malloc(lanesLen * sizeof(int32 *));
  int32  **oa = (int32 **)safe_malloc(lanesLen * sizeof(int32 *));

  for (int32 i=0; i<lanesLen; i++) {
    ia[i] = (int32 *)safe_calloc(cnsLen, sizeof(int32));
    oa[i] = (int32 *)safe_calloc(cnsLen, sizeof(int32));
  }


  for (int32 i=0; i<lanesLen; i++) {
    char *srow = multia[2*i];
    char *qrow = multia[2*i+1];

    assert(lanes[i].first != NULL);

    for (LaneNode *node=lanes[i].first; node != NULL; node = node->next) {
      int32 firstcol = (node->read->position.bgn < node->read->position.end) ? node->read->position.bgn : node->read->position.end;
      int32 lastcol  = (node->read->position.bgn < node->read->position.end) ? node->read->position.end : node->read->position.bgn;
      int32 orient   = (node->read->position.bgn < node->read->position.end) ? 1 : -1;

      //  Set ID and orientation

      for (int32 col=firstcol; col<lastcol; col++) {
        ia[i][col] = node->read->ident;
        oa[i][col] = orient;
      }

      //  Set bases

      int32 col  = firstcol;
      int32 cols = 0;

      for (int32 j=0; j<node->read->delta_length; j++) {
        int32 seglen = node->read->delta[j] - ((j > 0) ? node->read->delta[j-1] : 0);

        if (cols + seglen >= node->readLen)
          fprintf(stderr, "ERROR:  Clear ranges not correct.\n");
        assert(cols + seglen < node->readLen);

        memcpy(srow + col, node->sequence + cols, seglen);
        memcpy(qrow + col, node->quality  + cols, seglen);

        col += seglen;

        srow[col] = '-';
        qrow[col] = '-';
        col++;

        cols += seglen;
      }

      memcpy(srow + col, node->sequence + cols, node->readLen - cols);
      memcpy(qrow + col, node->quality  + cols, node->readLen - cols);
    }
  }

  //  Cleanup

  delete [] lanes;

  *array     = multia;
  *depth     = lanesLen;
  *id_array  = ia;
  *ori_array = oa;
}



void
PrintMultiAlignT(FILE *out,
                 MultiAlignT *ma,
                 gkStore *frag_store,
                 int32 show_qv,
                 int32 dots,
                 uint32 clrrng_flag)  {

  int32 depth;
  int32 i;
  int32 window;
  char **multia=NULL;
  int32 **idarray;
  int32 **oriarray;
  char *consensus = Getchar(ma->consensus,0);
  char *quality   = Getchar(ma->quality,0);

  gkFragment rsp;

  dots = 0;

  if ((consensus == NULL) || (consensus[0] == 0)) {
    fprintf(out, "No MultiAlignment to print for tig %d -- no consensus sequence present.\n", ma->maID);
    return;
  }

  int32 length = strlen(consensus);

  IMP2Array(GetIntMultiPos(ma->f_list,0),
            GetNumIntMultiPoss(ma->f_list),
            GetNumchars(ma->consensus),
            frag_store,
            &depth,
            &multia,
            &idarray,
            &oriarray,
            clrrng_flag);

  fprintf(out,"<<< begin Contig %d >>>",ma->maID);;

  char  gruler[MULTIALIGN_PRINT_WIDTH + 200];
  char  uruler[MULTIALIGN_PRINT_WIDTH + 200];

  int32 ungapped = 1;
  int32 tick     = 1;

  for (window=0;window<length;) {
    int32 row_id  = 0;
    int32 orient  = 0;
    int32 rowlen  = (window + MULTIALIGN_PRINT_WIDTH < length) ? MULTIALIGN_PRINT_WIDTH : length - window;

    fprintf(out, "\n");
    fprintf(out, "\n");
    fprintf(out, "<<<  Contig %d, gapped length: %d  >>>\n",ma->maID, length);

    {
      memset(gruler, 0, MULTIALIGN_PRINT_WIDTH + 200);
      memset(uruler, 0, MULTIALIGN_PRINT_WIDTH + 200);

      for (int32 rowind=0; rowind<rowlen; rowind++) {
        if (((window + 1 + rowind) % 25) == 0)
          sprintf(gruler + rowind, "| GAP=%d", window + 1 + rowind);

        if ((ungapped % 25) == 0)
          sprintf(uruler + rowind, "| UNG=%d", ungapped);

        if (consensus[window + rowind] != '-')
          ungapped++;
      }

      for (int32 i=0; i<MULTIALIGN_PRINT_WIDTH; i++) {
        if (gruler[i] == 0)
          gruler[i] = ' ';
        if (uruler[i] == 0)
          uruler[i] = ' ';
      }

      for (int32 i=MULTIALIGN_PRINT_WIDTH-1; (i >= 0) && (gruler[i] == ' '); i--)
        gruler[i] = 0;
      for (int32 i=MULTIALIGN_PRINT_WIDTH-1; (i >= 0) && (uruler[i] == ' '); i--)
        uruler[i] = 0;

      fprintf(out, "%s\n", gruler);
      fprintf(out, "%s\n", uruler);
    }


    {
      char save = consensus[window + rowlen];
      consensus[window+rowlen] = 0;
      fprintf(out,"%s  cns  (uid,iid) type\n", consensus+window);
      consensus[window+rowlen] = save;
    }

    {
      char save = quality[window + rowlen];
      quality[window+rowlen] = 0;
      fprintf(out,"%s  qlt\n", quality+window);
      quality[window+rowlen] = save;
    }

    for (i=0;i<depth;i++) {
      assert(multia[2*i] != NULL);

      //  Change matching bases to '.' or lowercase.
      //  Count the number of non-blank letters.

      int32  nonBlank = 0;

      for (int32 j=0; j<MULTIALIGN_PRINT_WIDTH; j++) {
        if (window + j > length)
          break;

        if (multia[2*i][window+j] == consensus[window+j]) {
          if (dots) {
            multia[2*i]  [window+j] = '.';
            multia[2*i+1][window+j] = ' ';
          } else {
            multia[2*i][window+j] = tolower(multia[2*i][window+j]);
          }
        }

        if (multia[2*i][window+j] != ' ')
          nonBlank++;

        if (idarray[i][window + j] > 0) {
          row_id = idarray[i][window + j];
          orient = oriarray[i][window + j];
        }
      }

      if (nonBlank == 0)
        continue;

      //  Figure out the ID and orientation for this block

      frag_store->gkStore_getFragment(row_id, &rsp, GKFRAGMENT_INF);

      {
        char save = multia[2*i][window + MULTIALIGN_PRINT_WIDTH];
        multia[2*i][window + MULTIALIGN_PRINT_WIDTH] = 0;
        fprintf(out, "%s   %c   (%s,%d)\n",
                multia[2*i]+window,
                (orient>0)?'>':'<',
                AS_UID_toString(rsp.gkFragment_getReadUID()),
                row_id);
        multia[2*i][window + MULTIALIGN_PRINT_WIDTH] = save;
      }

      if (show_qv) {
        char save = multia[2*i+1][window + MULTIALIGN_PRINT_WIDTH];
        multia[2*i+1][window + MULTIALIGN_PRINT_WIDTH] = 0;
        fprintf(out, "%s\n", multia[2*i+1]+window);
        multia[2*i+1][window + MULTIALIGN_PRINT_WIDTH] = save;
      }
    }

    window += MULTIALIGN_PRINT_WIDTH;
  }
  fprintf(out,"\n<<< end Contig %d >>>\n", ma->maID);



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
