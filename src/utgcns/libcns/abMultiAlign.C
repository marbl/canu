
#include "abAbacus.H"


void
abMultiAlign::getConsensus(abAbacus *abacus, char *bases, char *quals, uint32 &len, uint32 max) {
  abColumn   *col = abacus->getColumn(firstColumn());
  abBeadID   bid  = col->callID();

  len = 0;

  while (bid.isValid()) {
    abBead *bead = abacus->getBead(bid);

    bases[len] = abacus->getBase(bead->baseIdx());
    quals[len] = abacus->getQual(bead->baseIdx());

    //fprintf(stderr, "getConsensus()--  len %9u bead %5u base %c %c next %5d\n",
    //        len, bead->ident(), bases[len], quals[len], bead->nextID().get());

    bid = bead->nextID();
    len++;
  }

  bases[len] = 0;
  quals[len] = 0;

  if (len >= max) {
    display(abacus, stderr, 0, len+3);
    fprintf(stderr, "abMultiAlign::getConsensus()-- len=%u ! < max=%u\n", len, max);
  }
  assert(len < max);  //  Should be len+1 == max if we just reallocated, but max can be way bigger
}


void
abMultiAlign::getConsensus(abAbacus *abacus, tgTig *tig) {

  //fprintf(stderr, "abMultiAlign::getConsensus()-- transfering consensus of length %u to tig %u\n",
  //        length(), tig->tigID());

  //  +1 for the terminating nul byte.
  resizeArrayPair(tig->_gappedBases, tig->_gappedQuals, tig->_gappedLen, tig->_gappedMax, length() + 1, resizeArray_doNothing);

  //fprintf(stderr, "abMultiAlign::getConsensus()-- transfering consensus of length %u to tig %u (_gappedLen %u _gappedMax %u)\n",
  //        length(), tig->tigID(), tig->_gappedLen, tig->_gappedMax);

  getConsensus(abacus, tig->_gappedBases, tig->_gappedQuals, tig->_gappedLen, tig->_gappedMax);

  //fprintf(stderr, "abMultiAlign::getConsensus()-- transfering consensus of length %u to tig %u (_gappedLen %u _gappedMax %u)\n",
  //        length(), tig->tigID(), tig->_gappedLen, tig->_gappedMax);
}



static
uint32
getFragmentDeltas(abAbacus        *abacus,
                  abSeqID          sid, 
                  int32           *deltas) {

  abSequence         *seq = abacus->getSequence(sid);
  uint32              len = seq->length();

  abSeqBeadIterator  *fi  = abacus->createSeqBeadIterator(sid);
  uint32              dl  = 0;
  uint32              bp  = 0;

  //  index < length eliminates any endgaps from the delta list KAR, 09/19/02

  for (abBeadID bid=fi->next(); (bid.isValid()) && (bp < len); bid=fi->next()) {
    abBead  *bead = abacus->getBead(bid);

    if (abacus->getBase(bead->baseIdx()) == '-') {
      if (deltas)
        deltas[dl] = bp;
      dl++;
    }

    else
      bp++;
  }

  delete fi;

  if (deltas)
    deltas[dl] = 0;
  dl++;

  return(dl);
}



void
abMultiAlign::getPositions(abAbacus *abacus, tgTig *tig) {

  //  Interesting quandry.  We need to know the number of deltas so we can pre-allocate the array,
  //  or we need to keep growing it.  We're expecting lots of noise in the reads, and so lots of
  //  deltas, wich would mean lots of reallocating (and copying).  Is this faster than a pass
  //  through every read?

  uint32  nd = 0;

  for (abSeqID si=0; si.get()<abacus->numberOfSequences(); ++si)
    nd += getFragmentDeltas(abacus, si, NULL);

  resizeArray(tig->_childDeltas, tig->_childDeltasLen, tig->_childDeltasMax, nd, resizeArray_doNothing);

  tig->_childDeltasLen = 0;

  for (abSeqID si=0; si.get()<abacus->numberOfSequences(); ++si) {
    abSequence *seq = abacus->getSequence(si);

    //if (seq->multiAlignID() == abMultiAlignID())
    //  //  Not placed in the multialign??
    //  continue;
    //assert(seq->multiAlignID() == ident());

    int32 bgn = abacus->getColumn( (abacus->getBead( seq->firstBead() ))->colIdx())->position();
    int32 end = abacus->getColumn( (abacus->getBead( seq->lastBead()  ))->colIdx())->position();  // + seq->length() - 1))->colIdx())->position() + 1;

    tgPosition *child = tig->getChild(si.get());  //  Assumes one-to-one map of seqs to children;

    assert(seq->gkpIdent() == child->ident());

    if (seq->isRead() == true) {
      //child->anchor()      = 0;
      //child->aHang()       = 0;
      //child->bHang()       = 0;
      //child->bgnset()      = (seq->isForward() == true) ? bgn : end;
      //child->endset()      = (seq->isForward() == true) ? end : bgn;

      //child->_isReverse    = (seq->isForward() == false);

      //  The child remembers its orientation, and sets min/max appropriately.
      child->set(bgn, end);

      //  Grab the deltas for real now, and set the deltaLen to one less than the number returned
      //  (we don't care about the terminating zero!)

      child->_deltaOffset   = tig->_childDeltasLen;
      child->_deltaLen      = getFragmentDeltas(abacus, si, tig->_childDeltas + tig->_childDeltasLen);

      tig->_childDeltasLen += child->_deltaLen;

      child->_deltaLen--;

      //fprintf(stderr, "getPositions()-- si %u deltaLen %u offset %u\n",
      //        si.get(), child->_deltaLen, child->_deltaOffset);
    }
  }
}





void
abMultiAlign::display(abAbacus  *abacus,
                      FILE      *F,
                      uint32     from,
                      uint32     to) {

  if (to > length())
    to = length();

  if (from > to)
    return;

  int32   pageWidth = to - from;

  char   *sequence = new char [length() + 1];
  char   *quality  = new char [length() + 1];
  uint32  len      = 0;

  getConsensus(abacus, sequence, quality, len, length() + 1);

  if (len != length())
    fprintf(stderr, "abMultiAlign::display()-- len=%d != length=%d\n", len, length());
  //assert(len == length());


  uint32 numSeqs = abacus->numberOfSequences();

  abSeqBeadIterator     **fit       = new abSeqBeadIterator * [numSeqs];
  uint32                 *fid       = new uint32              [numSeqs];
  char                   *type      = new char                [numSeqs];  //  always 'r' for read
  uint32                 *bgn       = new uint32              [numSeqs];  //  former positions
  uint32                 *end       = new uint32              [numSeqs];



  for (uint32 i=0; i<numSeqs; i++) {
    abSequence  *seq = abacus->getSequence(i);

    //if (seq->multiAlignID() != ident()) {
    //  fprintf(stderr, "abMultiAlign::display()-- seq %d multialign %d != multialign %d\n", seq->ident().get(), seq->multiAlignID().get(), ident().get());
    //  continue;
    //}
    //assert(seq->multiAlignID() == ident());

    fid[i]  = seq->gkpIdent();
    fit[i]  = NULL;
    type[i] = (seq->isRead() == true) ? 'R' : '?';

    abColumn   *bgnColumn = abacus->getColumn(seq->firstBead());
    abColumn   *endColumn = abacus->getColumn(seq->lastBead());

    bgn[i] = bgnColumn->position();
    end[i] = endColumn->position();
  }

  uint32 window_start = from;

  fprintf(F,"\n==================== abMultiAlign::display %d ====================\n", ident().get());

  while (window_start < to) {
    fprintf(F,"\n");
    fprintf(F,"%d\n", window_start);
    fprintf(F,"%-*.*s <<< consensus\n", pageWidth, pageWidth, sequence + window_start);
    fprintf(F,"%-*.*s <<< quality\n\n", pageWidth, pageWidth, quality  + window_start);

    for (uint32 i=0; i<numSeqs; i++) {
      if (fid[i] == 0)
        continue;

      for (uint32 wi=window_start; wi<window_start + pageWidth; wi++) {


        //  If no valid iterator, 
        if (fit[i] == NULL) {

          if ((bgn[i] < wi) &&
              (wi     < end[i])) {
            //  Starting in the middle of the sequence, wi should equal window_start, right?

            //  This case only occurs if 'from' is set, and it's probably broken!

            //  It prints exactly one base, but serves to skip all the bases we shouldn't be printing.
            //  Once that base is printed, the iterator is valid, and we do the 'normal' case in the else clause.

            fit[i] = abacus->createSeqBeadIterator(i);

            abBeadID   bid = fit[i]->next();
            abColumn  *col = abacus->getColumn(bid);

            while (col->position() < wi) {
              bid = fit[i]->next();
              col = abacus->getColumn(bid);
            }

            fit[i]->prev();

          } else if (bgn[i] ==  wi) {
            //  Start at the beginning of the sequence
            fit[i] = abacus->createSeqBeadIterator(i);

          } else if ((window_start < bgn[i]) &&
                     (bgn[i]       < window_start+pageWidth)) {
            //  Spaces before the read
            fprintf(F," ");

          } else if ((window_start <= end[i]) &&
                     (end[i]       <  window_start+pageWidth)) {
            //  Spaces after the read
            fprintf(F," ");

          } else {
            //  Sequence isn't involved in this window.
            break;
          }
        }  //  End of fit[i] == NULL



        //  If a valid iterator, print a base.  If we've run out of bases, delete the iterator
        //  and print a space.
        //
        if (fit[i] != NULL) {
          abBeadID   bid = fit[i]->next();

          if (bid.isValid()) {
            char pc = abacus->getBase(bid);

            if (pc == sequence[wi])
              pc = tolower(pc);
            else
              pc = toupper(pc);

            fprintf(F, "%c", pc);
          }

          //  Not a valid bead, we just ran off the end of the read.  Print a space, and destroy the iterator.
          else {
            fprintf(F," ");
            delete fit[i];
            fit[i] = NULL;
          }
        }


        //  If at the end of the window, print the id of the object we just printed.
        //
        if (wi == window_start + pageWidth - 1)
          fprintf(F," <<< %d (%c)\n", fid[i], type[i]);
      }
    }

    window_start += pageWidth;
  }


  for (uint32 i=0; i<numSeqs; i++)
    delete fit[i];

  delete [] fit;
  delete [] fid;
  delete [] type;
  delete [] bgn;
  delete [] end;

  delete [] sequence;
  delete [] quality;
}
