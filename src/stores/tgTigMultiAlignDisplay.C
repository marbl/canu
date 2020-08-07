
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

#include "runtime.H"

#include "sqStore.H"
#include "tgStore.H"

#include "sequence.H"

#include <set>
#include <algorithm>

using namespace std;

//
//  The original utgcns consensus algorithm generated a gapped consensus sequence, which
//  made displaying a multialignment of reads The original version only worked with a gapped consensus sequence.
//

class alignRowEntry {
public:
  alignRowEntry(sqStore *seq_, tgTig *tig_, uint32 child_);

  ~alignRowEntry() {
    delete [] bases;
    delete [] quals;
    delete [] delta;
  };

  tgPosition     *position;
  int32           sequenceLength;

  char           *bases;
  char           *quals;

  int32          *delta;
  int32           deltaLength;

  alignRowEntry  *nextEntry;
};



alignRowEntry::alignRowEntry(sqStore     *seq_,
                             tgTig       *tig_,
                             uint32       child_) {

  sqRead     *read = seq_->sqStore_getRead(tig_->getChild(child_)->ident(), new sqRead());

  //  Set basic stuff and allocate space.

  position        = tig_->getChild(child_);
  sequenceLength  = read->sqRead_length();

  bases           = new char  [sequenceLength + 1];
  quals           = new char  [sequenceLength + 1];
  delta           = new int32 [position->deltaLength()];
  deltaLength     = position->deltaLength();

  nextEntry       = NULL;

  //  Copy sequence and create empty quals.

  char  *b = read->sqRead_sequence();

  for (uint32 ii=0; ii<sequenceLength; ii++) {
    bases[ii] = b[ii];
    quals[ii] = '!';
  }

  bases[sequenceLength] = 0;
  quals[sequenceLength] = 0;

  //  Flip the sequence if needed.  The deltas are relative to this flipped sequence.

  if (position->isReverse())
    ::reverseComplement(bases, quals, sequenceLength);

  //  Decode deltas.

  tig_->_childDeltaBits->setPosition(position->deltaOffset());

  for (uint32 dd=0; dd<deltaLength; dd++) {
    delta[dd] = tig_->_childDeltaBits->getEliasDelta();

    if (tig_->_childDeltaBits->getBit() == 0)
      delta[dd] = -delta[dd];
  }

  //  Cleanup.

  delete read;
}



class alignRow {
public:
  alignRow() {
    firstEntry = NULL;
    lastEntry  = NULL;
    lastColumn = 0;
  };

  ~alignRow() {
    alignRowEntry  *e = firstEntry;
    alignRowEntry  *n = NULL;

    while (e) {
      n = e->nextEntry;
      delete e;
      e = n;
    }
  };

  bool    addEntry(alignRowEntry *entry, uint32 spacing);

  alignRowEntry   *firstEntry;
  alignRowEntry   *lastEntry;

  int32            lastColumn;
};



bool
alignRow::addEntry(alignRowEntry *entry, uint32 spacing) {
  int32 leftPos = entry->position->min();

  if ((lastColumn > 0) &&                  //  If something in this row, fail
      (leftPos < lastColumn + spacing))    //  if the new entry intersects.
    return(false);

  assert(entry->nextEntry == NULL);

  if (firstEntry == NULL) {                //  Add the entry to our list
    firstEntry           = entry;          //  of entries.
    lastEntry            = entry;
  } else {
    lastEntry->nextEntry = entry;
    lastEntry            = entry;
    lastEntry->nextEntry = NULL;
  }

  lastColumn = leftPos + entry->deltaLength + entry->sequenceLength;

  return(true);
}







void
tgTig::display(FILE     *F,
               sqStore  *seq,
               uint32    displayWidth,
               uint32    displaySpacing,
               bool      withQV,
               bool      withDots)  {

  withQV   = false;
  withDots = true;

  if (consensusExists() == 0) {
    fprintf(F, "No MultiAlignment to print for tig %d -- no consensus sequence present.\n", tigID());
    return;
  }

  fprintf(stderr, "tgTig::display()--  display tig %d with %d children\n", tigID(), _childrenLen);
  fprintf(stderr, "tgTig::display()--  width %u spacing %u\n", displayWidth, displaySpacing);

  //
  //  Convert the children to a list of rows to print.  Worst case is one row per child.
  //

  alignRow  *row     = NULL;
  alignRow  *rows    = new alignRow [_childrenLen];
  int32      rowsLen = 0;

  // Sort the children by leftmost position within tig.

  std::sort(_children, _children + _childrenLen);

  //  Load into rows.

  for (int32 i=0; i<_childrenLen; i++) {
    alignRowEntry  *entry = new alignRowEntry(seq, this, i);
    uint32          rowsPos;

    //  Try to add this new node to the rows.  The last iteration will always succeed, adding the
    //  node to a fresh empty lane.

    for (rowsPos=0; rowsPos <= rowsLen; rowsPos++)
      if (rows[rowsPos].addEntry(entry, displaySpacing))
        break;

    assert(rowsPos <= rowsLen);

    //  If it is added to the last one, increment our cnsLen.

    if (rowsPos == rowsLen)
      rowsLen++;
  }

  //  Find the list of gaps we need to insert into the reference.

  set<int32>   gapPositions;

#if 0
  //  This is wrong, each delta is relative to the last one.
  for (int32 rr=0; rr<rowsLen; rr++)
    for (alignRowEntry *entry=rows[rr].firstEntry; entry != NULL; entry = entry->nextEntry)
      for (uint32 dd=0; dd<entry->deltaLength; dd++)
        if (entry->delta[dd] < 0)
          gapPositions.insert(entry->position->min() + -entry->delta[dd] - 1);
#endif

  //  Allocate space.

  char   **displayBases  = new char  * [rowsLen];   //  Bases to print in the current window.
  char   **displayQuals  = new char  * [rowsLen];   //  Quals to print in the current window.
  int32  **displayIDs    = new int32 * [rowsLen];   //  The ID present at this position.
  char   **displayFwd    = new char  * [rowsLen];   //

  for (uint32 ii=0; ii<rowsLen; ii++) {
    displayBases[ii] = new char  [length() + gapPositions.size() + displayWidth + 1];
    displayQuals[ii] = new char  [length() + gapPositions.size() + displayWidth + 1];
    displayIDs[ii]   = new int32 [length() + gapPositions.size() + 1];
    displayFwd[ii]   = new char  [length() + gapPositions.size() + 1];

    memset(displayBases[ii], ' ', sizeof(char)  * (length() + gapPositions.size() + displayWidth + 1));
    memset(displayQuals[ii], ' ', sizeof(char)  * (length() + gapPositions.size() + displayWidth + 1));

    memset(displayIDs[ii],    0,  sizeof(int32) * (length() + gapPositions.size() + 1));
    memset(displayFwd[ii],    0,  sizeof(char)  * (length() + gapPositions.size() + 1));
  }

  //  Copy read sequence to the display arrays, ignoring gaps in the tig for now.
  //
  //  Set bases.  Three cases.
  //    No gap      - just copy the base.
  //    Gap in read - set to 'gap' base/qual
  //    Gap in tig  - ignore for now

  for (int32 rr=0; rr<rowsLen; rr++) {
    for (alignRowEntry *entry=rows[rr].firstEntry; entry != NULL; entry = entry->nextEntry) {
      int32 readPos   = 0;
      int32 activeCol = entry->position->min();

      for (int32 j=0; j<entry->deltaLength; j++) {
        int32 delta  = entry->delta[j];
        int32 seglen = ((delta > 0) ? delta : -delta) - 1;

        //  Copy bases from the read to the alignment.  We ignore gaps in the tig caused
        //  by other reads (gapPositions.count(activeCol+cc) > 0) here.

        for (int32 cc=0; cc<seglen; cc++) {
          displayBases[rr][activeCol]   = entry->bases[readPos];
          displayQuals[rr][activeCol]   = entry->quals[readPos++];
          displayIDs  [rr][activeCol]   = entry->position->ident();
          displayFwd  [rr][activeCol++] = entry->position->isForward();
        }

        //  Add a gap into the read.

        if (delta > 0) {
          displayBases[rr][activeCol]   = '-';
          displayQuals[rr][activeCol]   = '-';
          displayIDs  [rr][activeCol]   = entry->position->ident();
          displayFwd  [rr][activeCol++] = entry->position->isForward();
        }

        //  But ignore gaps into the tig, and skip over the base in the read.
        else {
          readPos++;
        }
      }

      //  Copy in the rest of the sequence.
#if 1
      while (readPos < entry->sequenceLength) {
        displayBases[rr][activeCol]   = entry->bases[readPos];
        displayQuals[rr][activeCol]   = entry->quals[readPos++];
        displayIDs  [rr][activeCol]   = entry->position->ident();
        displayFwd  [rr][activeCol++] = entry->position->isForward();
      }
#endif
      //memcpy(srow + col, entry->bases + cols, entry->readLen - cols);
      //memcpy(qrow + col, entry->quals + cols, entry->readLen - cols);
    }
  }

  //  Cleanup

  delete [] rows;


  //
  //
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


    {
      char save = _bases[window + rowlen];
      _bases[window + rowlen] = 0;

      fprintf(F, "%s  cns  (iid) type\n", _bases + window);

      _bases[window + rowlen] = save;
    }

    {
      for (uint32 ii=0; ii<rowlen; ii++)   //  Adjust QV for display.
        _quals[window+ii] += '!';

      char save = _quals[window + rowlen];
      _quals[window+rowlen] = 0;     //  Terminate the substring.

      fprintf(F, "%s  qlt\n", _quals + window);

      _quals[window+rowlen] = save;  //  Unterminate.

      for (uint32 ii=0; ii<rowlen; ii++)   //  Unadjust.
        _quals[window+ii] -= '!';
    }

    //  Display.

    for (uint32 i=0; i<rowsLen; i++) {
      int32 row_id = -1;
      bool  isfwd = false;

      //  Change matching bases to '.' or lowercase.
      //  Count the number of non-blank letters.

      for (int32 j=0; j<displayWidth; j++) {
        if (window + j > length())
          break;

        if (displayBases[i][window+j] == _bases[window+j]) {
          if (withDots) {
            displayBases[i][window+j] = '.';
            displayQuals[i][window+j] = ' ';
          } else {
            displayBases[i][window+j] = tolower(displayBases[i][window+j]);
          }
        }

        if (displayIDs[i][window + j] > 0) {
          row_id = displayIDs[i][window + j];
          isfwd  = displayFwd[i][window + j];
        }
      }

      if (row_id == -1)
        continue;

      //  Display the bases in this row, with orientation and the id.

      {
        char save = displayBases[i][window + displayWidth];
        displayBases[i][window + displayWidth] = 0;

        fprintf(F, "%s   %c   (%d)\n", displayBases[i] + window, (isfwd) ? '>' : '<', row_id);

        displayBases[i][window + displayWidth] = save;
      }

      //  Display the quals in this row.

      if (withQV) {
        char save = displayQuals[i][window + displayWidth];
        displayQuals[i][window + displayWidth] = 0;

        fprintf(F, "%s\n", displayQuals[i] + window);

        displayQuals[i][window + displayWidth] = save;
      }
    }

    window += displayWidth;
  }

  fprintf(F, "\n<<< end Contig %d >>>\n", tigID());

  delete [] uruler;
  delete [] gruler;

  for (uint32 i=0; i < rowsLen; i++) {
    delete [] displayBases[i];
    delete [] displayQuals[i];
    delete [] displayIDs[i];
    delete [] displayFwd[i];
  }

  delete [] displayBases;
  delete [] displayQuals;
  delete [] displayIDs;
  delete [] displayFwd;

  //fprintf(stderr, "Return.\n");
}
