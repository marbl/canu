
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


#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"

#include "AS_UTL_reverseComplement.H"

#include <algorithm>


class LaneNode {
public:
  LaneNode() {
    read      = NULL;
    readLen   = 0;
    bases     = NULL;
    quals     = NULL;
    next      = NULL;
  };

  ~LaneNode() {
    delete [] bases;
    delete [] quals;
  };

  tgPosition       *read;
  int32             readLen;

  char             *bases;   //  Allocated bases
  char             *quals;   //  Allocated quals
  int32            *delta;   //  Pointer to tgTig _childDeltas

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
  bool addNode(LaneNode *node, uint32 displaySpacing) {
    int32 leftpos = node->read->min();

    if ((lastcol > 0) &&
        (leftpos < lastcol + displaySpacing))
      return(false);

    assert(node->next == NULL);

    if (first == NULL) {
      first      = node;
      last       = node;
    } else {
      last->next = node;
      last       = node;
    }

    lastcol = leftpos + node->read->deltaLength() + node->readLen;

    return(true);
  };

  LaneNode         *first;
  LaneNode         *last;
  int32             lastcol;
};





void
tgTig::display(FILE     *F,
               gkStore  *gkp,
               uint32    displayWidth,
               uint32    displaySpacing,
               bool      withQV,
               bool      withDots)  {

  uint32     depth;
  char     **multia=NULL;
  int32    **idarray;
  int32    **oriarray;

  if (gappedLength() == 0) {
    fprintf(F, "No MultiAlignment to print for tig %d -- no consensus sequence present.\n", tigID());
    return;
  }

  //
  //  Convert the children to a list of lines to print
  //

  Lane  *lane     = NULL;
  Lane  *lanes    = new Lane [_childrenLen];
  int32  lanesLen = 0;
  int32  lanesPos = 0;

  // Sort the fragments by leftmost position within tig

  std::sort(_children, _children + _childrenLen);

  //  Load into lanes.

  for (int32 i=0; i<_childrenLen; i++) {
    gkRead     *read      = gkp->gkStore_getRead(_children[i].ident());  //  Too many reads in this code.
    gkReadData  readData;

    gkp->gkStore_loadReadData(read, &readData);

    LaneNode   *node = new LaneNode();

    node->read      = _children + i;
    node->readLen   = read->gkRead_clearRegionLength();

    node->bases     = new char [node->readLen + 1];
    node->quals     = new char [node->readLen + 1];

    node->delta     = _childDeltas + node->read->deltaOffset();

    memcpy(node->bases, readData.gkReadData_getSequence()  + read->gkRead_clearRegionBegin(), sizeof(char) * node->readLen);
    memcpy(node->quals, readData.gkReadData_getQualities() + read->gkRead_clearRegionBegin(), sizeof(char) * node->readLen);

    node->bases[node->readLen] = 0;
    node->quals [node->readLen] = 0;

    if (node->read->isReverse())
      reverseComplement(node->bases, node->quals, node->readLen);

    //  Try to add this new node to the lanes.  The last iteration will always succeed, adding the
    //  node to a fresh empty lane.

    for (lanesPos=0; lanesPos <= lanesLen; lanesPos++)
      if (lanes[lanesPos].addNode(node, displaySpacing))
        break;

    assert(lanesPos <= lanesLen);

    //  If it is added to the last one, increment our cnsLen.

    if (lanesPos == lanesLen)
      lanesLen++;
  }

  //  Process.

  multia = new char * [2 * lanesLen];

  for (int32 i=0; i<2*lanesLen; i++) {
    multia[i]  = new char [gappedLength() + 1 + displayWidth];

    memset(multia[i], ' ', gappedLength() + displayWidth);
    multia[i][gappedLength()] = 0;
  }

  idarray  = new int32 * [lanesLen];
  oriarray = new int32 * [lanesLen];

  for (int32 i=0; i<lanesLen; i++) {
    idarray[i] = new int32 [gappedLength()];
    oriarray[i] = new int32 [gappedLength()];

    memset(idarray[i], 0, sizeof(int32) * gappedLength());
    memset(oriarray[i], 0, sizeof(int32) * gappedLength());
  }


  for (int32 i=0; i<lanesLen; i++) {
    char *srow = multia[2*i];
    char *qrow = multia[2*i+1];

    assert(lanes[i].first != NULL);

    for (LaneNode *node=lanes[i].first; node != NULL; node = node->next) {
      int32 firstcol = (node->read->bgn() < node->read->end()) ? node->read->bgn() : node->read->end();
      int32 lastcol  = (node->read->bgn() < node->read->end()) ? node->read->end() : node->read->bgn();
      int32 orient   = (node->read->bgn() < node->read->end()) ? 1 : -1;

      //  Set ID and orientation

      for (int32 col=firstcol; col<lastcol; col++) {
        idarray[i][col] = node->read->ident();
        oriarray[i][col] = orient;
      }

      //  Set bases

      int32 col  = firstcol;
      int32 cols = 0;

      for (int32 j=0; j<node->read->deltaLength(); j++) {
        int32 seglen = node->delta[j] - ((j > 0) ? node->delta[j-1] : 0);

        if (cols + seglen >= node->readLen)
          fprintf(stderr, "ERROR:  Clear ranges not correct.\n");
        assert(cols + seglen < node->readLen);

        memcpy(srow + col, node->bases + cols, seglen);
        memcpy(qrow + col, node->quals  + cols, seglen);

        col += seglen;

        srow[col] = '-';
        qrow[col] = '-';
        col++;

        cols += seglen;
      }

      memcpy(srow + col, node->bases + cols, node->readLen - cols);
      memcpy(qrow + col, node->quals  + cols, node->readLen - cols);
    }
  }

  //  Cleanup

  delete [] lanes;


  //
  //
  //

  fprintf(F, "<<< begin Contig %d >>>", tigID());

  char *gruler = new char [displayWidth + 200];
  char *uruler = new char [displayWidth + 200];

  int32 ungapped = 1;
  int32 tick     = 1;

  for (uint32 window=0; window < gappedLength(); ) {
    uint32 row_id  = 0;
    uint32 orient  = 0;
    uint32 rowlen  = (window + displayWidth < gappedLength()) ? displayWidth : gappedLength() - window;

    fprintf(F, "\n");
    fprintf(F, "\n");
    fprintf(F, "<<<  tig %d, gapped length: %d  >>>\n", tigID(), gappedLength());

    {
      memset(gruler, 0, displayWidth + 200);
      memset(uruler, 0, displayWidth + 200);

      for (uint32 rowind=0; rowind<rowlen; rowind++) {
        if (((window + 1 + rowind) % 25) == 0)
          sprintf(gruler + rowind, "| GAP=%d", window + 1 + rowind);

        if ((ungapped % 25) == 0)
          sprintf(uruler + rowind, "| UNG=%d", ungapped);

        if (_gappedBases[window + rowind] != '-')
          ungapped++;
      }

      for (int32 i=0; i<displayWidth; i++) {
        if (gruler[i] == 0)
          gruler[i] = ' ';
        if (uruler[i] == 0)
          uruler[i] = ' ';
      }

      for (int32 i=displayWidth-1; (i >= 0) && (gruler[i] == ' '); i--)
        gruler[i] = 0;
      for (int32 i=displayWidth-1; (i >= 0) && (uruler[i] == ' '); i--)
        uruler[i] = 0;

      fprintf(F, "%s\n", gruler);
      fprintf(F, "%s\n", uruler);
    }


    {
      char save = _gappedBases[window + rowlen];
      _gappedBases[window + rowlen] = 0;

      fprintf(F, "%s  cns  (iid) type\n", _gappedBases + window);

      _gappedBases[window + rowlen] = save;
    }

    {
      char save = _gappedQuals[window + rowlen];
      _gappedQuals[window+rowlen] = 0;

      fprintf(F, "%s  qlt\n", _gappedQuals + window);

      _gappedQuals[window+rowlen] = save;
    }

    for (uint32 i=0; i<depth; i++) {
      assert(multia[2*i] != NULL);

      //  Change matching bases to '.' or lowercase.
      //  Count the number of non-blank letters.

      int32  nonBlank = 0;

      for (int32 j=0; j<displayWidth; j++) {
        if (window + j > gappedLength())
          break;

        if (multia[2*i][window+j] == _gappedBases[window+j]) {
          if (withDots) {
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

      {
        char save = multia[2*i][window + displayWidth];
        multia[2*i][window + displayWidth] = 0;

        fprintf(F, "%s   %c   (%d)\n", multia[2*i] + window, (orient>0)?'>':'<', row_id);

        multia[2*i][window + displayWidth] = save;
      }

      if (withQV) {
        char save = multia[2*i+1][window + displayWidth];
        multia[2*i+1][window + displayWidth] = 0;

        fprintf(F, "%s\n", multia[2*i+1] + window);

        multia[2*i+1][window + displayWidth] = save;
      }
    }

    window += displayWidth;
  }

  fprintf(F, "\n<<< end Contig %d >>>\n", tigID());

  for (uint32 i=0; i < 2*depth; i++)
    delete [] multia[i];

  delete [] multia;

  for (uint32 i=0; i < depth; i++) {
    delete [] idarray[i];
    delete [] oriarray[i];
  }

  delete [] idarray;
  delete [] oriarray;
}