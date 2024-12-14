
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */


#include "sqStore.H"
#include "tgStore.H"

#include "sequence.H"

#include <set>
#include <algorithm>

//
//  The original utgcns consensus algorithm generated a gapped consensus sequence, which
//  made displaying a multialignment of reads The original version only worked with a gapped consensus sequence.
//

class alignRowEntry {
public:
  alignRowEntry(sqStore *seq_,
                tgTig   *tig_,
                uint32   child_) {
    sqRead  read;

    seq_->sqStore_getRead(tig_->getChild(child_)->ident(), &read);

    position        = tig_->getChild(child_);          //  Set basic stuff and allocate space.
    sequenceLength  = read.sqRead_length();

    bases           = new char  [sequenceLength + 1];    bases[sequenceLength] = 0;
    quals           = new char  [sequenceLength + 1];    quals[sequenceLength] = 0;
    delta           = new int32 [position->deltaLength()];
    deltaLength     =            position->deltaLength();

    char  *b = read.sqRead_sequence();                 //  Copy sequence and create empty quals.

    for (uint32 ii=0; ii<sequenceLength; ii++) {
      bases[ii] = b[ii];
      quals[ii] = '!';
    }

    if (position->isReverse())                         //  Reverse complement if needed.
      ::reverseComplement(bases, quals, sequenceLength);

    tig_->_childDeltaBits->setPosition(position->deltaOffset());

    for (uint32 dd=0; dd<deltaLength; dd++) {          //  Decode deltas.
      delta[dd] = tig_->_childDeltaBits->getEliasDelta();

      if (tig_->_childDeltaBits->getBit() == 0)
        delta[dd] = -delta[dd];

      //fprintf(stderr, "READ %8u DELTA[%02u] %d\n", tig_->getChild(child_)->ident(), dd, delta[dd]);
    }
  }

  ~alignRowEntry() {
    delete [] bases;
    delete [] quals;
    delete [] delta;
  }

  tgPosition     *position       = nullptr;
  int32           sequenceLength = 0;

  char           *bases          = nullptr;
  char           *quals          = nullptr;

  int32          *delta          = nullptr;
  int32           deltaLength    = 0;

  alignRowEntry  *nextEntry      = nullptr;
};


class alignRow {
public:
  alignRow() {
  }

  ~alignRow() {
    for (alignRowEntry *e=firstEntry, *n; e; e=n) {   //  Walk down the singly-linked list,
      n = e->nextEntry;                               //  deleting nodes and moving to the next.
      delete e;
    }
  }

  bool    addEntry(alignRowEntry *entry, uint32 spacing) {
    int32 leftPos = entry->position->min();

    if ((lastColumn > 0) &&                  //  If something in this row, fail
        (leftPos < lastColumn + spacing))    //  if the new entry intersects.
      return false;

    assert(entry->nextEntry == nullptr);

    if (firstEntry == nullptr) {             //  Add the entry to our list
      firstEntry           = entry;          //  of entries.
      lastEntry            = entry;
    } else {
      lastEntry->nextEntry = entry;
      lastEntry            = entry;
    }

    lastColumn = leftPos + entry->deltaLength + entry->sequenceLength;

    return true;
  }

  alignRowEntry   *firstEntry = nullptr;
  alignRowEntry   *lastEntry  = nullptr;

  int32            lastColumn = 0;
};





void
tgTig::display(FILE     *F,
               sqStore  *seq,
               uint32    displayWidth,
               uint32    displaySpacing,
               bool      withDots)  {

  if (consensusExists() == false) {
    fprintf(F, "No MultiAlignment to print for tig %d -- no consensus sequence present.\n", tigID());
    return;
  }

  // Sort the children by leftmost position within tig.

  std::sort(_children, _children + _childrenLen);

  //  Assign children to rows, reusing rows when they have empty space.
  //    row1 -- [read1  read4  .]
  //    row2 -- [   read2  .....]
  //    row3 -- [     read3  ...]
  //  Worst case is one row per read.
  //
  //  The inner loop tries to add the read ('entry') to each existing row.
  //  If there is space for the read, true is returned and the loop exits.
  //  If all existing rows are full, a fresh empty row is used, and the
  //  number of active rows is increased.

  alignRow  *row     = nullptr;
  alignRow  *rows    = new alignRow [_childrenLen];
  uint32     rowsLen = 0;

  for (int32 i=0; i<_childrenLen; i++) {
    alignRowEntry  *entry = new alignRowEntry(seq, this, i);   //  Everything we need to know about a read in a row.

    for (uint32 rr=0; rr<_childrenLen; rr++)
      if (rows[rr].addEntry(entry, displaySpacing)) {
        rowsLen = std::max(rowsLen, rr+1);
        break;
      }
  }

  //  Find the list of gaps we need to insert into the reference.

  uint32  gapPositionsSize = 1024;

#if 0
  std::set<int32>   gapPositions;

  //  This is wrong, each delta is relative to the last one.
  for (int32 rr=0; rr<rowsLen; rr++)
    for (alignRowEntry *entry=rows[rr].firstEntry; entry != nullptr; entry = entry->nextEntry)
      for (uint32 dd=0; dd<entry->deltaLength; dd++)
        if (entry->delta[dd] < 0)
          gapPositions.insert(entry->position->min() + -entry->delta[dd] - 1);
#endif

  //  Allocate space for the display - a matrix of [rows][columns], where the
  //  number of rows is as above, and the number of columns is the tigLength
  //  + gaps + a little bit.

  char   **displayBases  = new char  * [rowsLen];   //  Base at this [row][column]
  int32  **displayIDs    = new int32 * [rowsLen];   //  Read ID
  char   **displayFwd    = new char  * [rowsLen];   //  Read orientation

  for (uint32 rr=0; rr<rowsLen; rr++) {
    uint32  sl = length() + gapPositionsSize + displayWidth + 1;   //  Sequence length
    uint32  il = length() + gapPositionsSize + displayWidth + 1;   //  Information length

    displayBases[rr] = new char  [sl];    memset(displayBases[rr], ' ', sizeof(char)  * sl);
    displayIDs[rr]   = new int32 [il];    memset(displayIDs[rr],    0,  sizeof(int32) * il);
    displayFwd[rr]   = new char  [il];    memset(displayFwd[rr],    0,  sizeof(char)  * il);
  }

  //  Copy read data to the display arrays.  We cannot show gaps in the tig.
  //
  //  Each row has a singly-linked list of the reads that are contained it.

  for (int32 rr=0; rr<rowsLen; rr++) {
    for (alignRowEntry *entry=rows[rr].firstEntry; entry != nullptr; entry = entry->nextEntry) {
      int32 rPos = 0;                        //  Position we're at in the read.
      int32 dPos = entry->position->min();   //  Position we're at in the row.

      for (int32 j=0; j<entry->deltaLength; j++) {
        int32 delta  = entry->delta[j];
        int32 seglen = ((delta > 0) ? delta : -delta) - 1;

        //  Copy seglen bases from the read to the display.

        for (int32 cc=0; cc<seglen; cc++, rPos++, dPos++) {
          displayBases[rr][dPos] = entry->bases[rPos];
          displayIDs  [rr][dPos] = entry->position->ident();
          displayFwd  [rr][dPos] = entry->position->isForward();
        }

        //  Add a gap.
        //    delta > 0 -> gap in the read, extra base in the tig
        //    delta < 0 -> gap in the tig, extra base in the read - these are ignored

        if (delta > 0) {
          displayBases[rr][dPos] = '-';
          displayIDs  [rr][dPos] = entry->position->ident();
          displayFwd  [rr][dPos] = entry->position->isForward();
          dPos++;
        }
        else {
          rPos++;
        }
      }

      //  Done with the gaps.  Copy the rest of the sequence.

      for (int32 cc=0; rPos < entry->sequenceLength; rPos++, dPos++) {
        displayBases[rr][dPos] = entry->bases[rPos];
        displayIDs  [rr][dPos] = entry->position->ident();
        displayFwd  [rr][dPos] = entry->position->isForward();
      }
    }
  }

  //
  //  DISPLAY!
  //

  fprintf(F, "<<< begin Contig %d >>>", tigID());

  uint32  lruler = displayWidth + 200;
  char   *gruler = new char [lruler];
  char   *uruler = new char [lruler];

  int32 ungapped = 1;
  int32 tick     = 1;

  //  The original used 'length = strlen(consensus)' which does NOT include the terminating NUL.

  for (uint32 window=0; window < length(); ) {
    uint32 rowlen  = (window + displayWidth < length()) ? displayWidth : length() - window;

    fprintf(F, "\n");
    fprintf(F, "\n");
    fprintf(F, "<<<  tig %d, gapped length: %d  >>>\n", tigID(), length());

    {
      memset(gruler, 0, displayWidth + 200);
      memset(uruler, 0, displayWidth + 200);

      for (uint32 rowind=0; rowind<rowlen; rowind++) {
        if (((window + 1 + rowind) % 25) == 0)
          snprintf(gruler + rowind, lruler, "| GAP=%d", window + 1 + rowind);

        if ((ungapped % 25) == 0)
          snprintf(uruler + rowind, lruler, "| UNG=%d", ungapped);

        if (_bases[window + rowind] != '-')
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


    if (1) {  //  Display the tig sequence.
      char save = _bases[window + rowlen];
      _bases[window + rowlen] = 0;

      fprintf(F, "%s  cns  (iid) type\n", _bases + window);

      _bases[window + rowlen] = save;
    }

    if (0) {  //  Display the tig quality.
      for (uint32 ii=0; ii<rowlen; ii++)   //  Adjust QV for display.
        _quals[window+ii] += '!';

      char save = _quals[window + rowlen];
      _quals[window+rowlen] = 0;     //  Terminate the substring.

      fprintf(F, "%s  qlt\n", _quals + window);

      _quals[window+rowlen] = save;  //  Unterminate.

      for (uint32 ii=0; ii<rowlen; ii++)   //  Unadjust.
        _quals[window+ii] -= '!';
    }

    //  Display each row.

    for (uint32 rr=0; rr<rowsLen; rr++) {
      int32 row_id = -1;
      bool  isfwd = false;

      //  Change matching bases to '.' or lowercase.

      for (int32 cc=0; cc<displayWidth; cc++) {
        if (window + cc > length())
          break;

        if (displayBases[rr][window+cc] == _bases[window+cc]) {
          if (withDots) {
            displayBases[rr][window+cc] = '.';
          } else {
            displayBases[rr][window+cc] = tolower(displayBases[rr][window+cc]);
          }
        }

        if (displayIDs[rr][window + cc] > 0) {
          row_id = displayIDs[rr][window + cc];
          isfwd  = displayFwd[rr][window + cc];
        }
      }

      if (row_id == -1)  //  Nothing in this row
        continue;

      //  Display the bases in this row, with orientation and the id.

      {
        char save = displayBases[rr][window + displayWidth];
        displayBases[rr][window + displayWidth] = 0;

        fprintf(F, "%s   %c   (%d)\n", displayBases[rr] + window, (isfwd) ? '>' : '<', row_id);

        displayBases[rr][window + displayWidth] = save;
      }
    }

    window += displayWidth;
  }

  fprintf(F, "\n<<< end Contig %d >>>\n", tigID());

  delete [] uruler;
  delete [] gruler;

  for (uint32 rr=0; rr<rowsLen; rr++) {
    delete [] displayBases[rr];
    delete [] displayIDs[rr];
    delete [] displayFwd[rr];
  }

  delete [] displayBases;
  delete [] displayIDs;
  delete [] displayFwd;

  //fprintf(stderr, "Return.\n");
}
