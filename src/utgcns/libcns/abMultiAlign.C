
#include "abAbacus.H"


void
abMultiAlign::getConsensus(abAbacus *abacus, char *bases, char *quals, uint32 &len, uint32 max) {
  abColumn   *col = abacus->getColumn(firstColumn());
  abBeadID   bid  = col->callID();

  for (len=0; bid.isValid(); len++) {
    abBead *bead = abacus->getBead(bid);

    bases[len] = abacus->getBase(bead->baseIdx());
    quals[len] = abacus->getQual(bead->baseIdx());

    //fprintf(stderr, "getConsensus()--  len %9u bead %5u base %c %c next %5d\n",
    //        len, bead->ident(), bases[len], quals[len], bead->nextID().get());

    bid = bead->nextID();
  }

  bases[len] = 0;
  quals[len] = 0;

  assert(len < max);  //  Should be len+1 == max if we just reallocated, but max can be way bigger
}


void
abMultiAlign::getConsensus(abAbacus *abacus, tgTig *tig) {

  fprintf(stderr, "abMultiAlign::getConsensus()-- transfering consensus of length %u to tig %u\n", length(), tig->tigID());

  //  +1 for the terminating nul byte.
  resizeArrayPair(tig->_gappedBases, tig->_gappedQuals, tig->_gappedLen, tig->_gappedMax, length() + 1, resizeArray_doNothing);

  getConsensus(abacus, tig->_gappedBases, tig->_gappedQuals, tig->_gappedLen, tig->_gappedMax);
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

    if (abacus->getBase(bead->baseIdx()) == '-')
      if (deltas)
        deltas[dl++] = bp;

    else
      bp++;
  }

  return(dl);
}



void
abMultiAlign::getPositions(abAbacus *abacus, tgTig *tig) {


  //  Interesting quandry.  We need to know the number of deltas so we can pre-allocate the array,
  //  or we need to keep growing it.  We're expecting lots of noise in the reads, and so lots of
  //  deltas, wich would mean lots of reallocating (and copying).  Is this faster than a pass
  //  through every read?

  uint32  nd = 0;

  for (abSeqID si=0; si.get()<abacus->numberOfSequences(); ++si) {
    nd += getFragmentDeltas(abacus, si, NULL);
  }

  resizeArray(tig->_childDeltas, tig->_childDeltasLen, tig->_childDeltasMax, nd, resizeArray_doNothing);

  tig->_childDeltasLen = 0;


  for (abSeqID si=0; si.get()<abacus->numberOfSequences(); ++si) {
    abSequence *seq = abacus->getSequence(si);

    if (seq->multiAlignID() == abMultiAlignID())
      //  Not placed in the multialign??
      continue;

    assert(seq->multiAlignID() == ident());

    uint32 bgn = abacus->getColumn( (abacus->getBead( seq->firstBead()                    ))->colIdx())->position();
    uint32 end = abacus->getColumn( (abacus->getBead( seq->firstBead() + seq->length() - 1))->colIdx())->position() + 1;

    tgPosition *child = tig->getChild(si.get());  //  Assumes one-to-one map of seqs to children;

    assert(seq->gkpIdent() == child->ident());

    if (seq->isRead() == true) {
      child->anchor()      = 0;
      child->aHang()       = 0;
      child->bHang()       = 0;
      child->bgn()         = (seq->isForward() == true) ? bgn : end;
      child->end()         = (seq->isForward() == true) ? end : bgn;

      child->_deltaOffset  = tig->_childDeltasLen;
      child->_deltaLen     = getFragmentDeltas(abacus, si, tig->_childDeltas + tig->_childDeltasLen);;

      tig->_childDeltasLen += tig->_childDeltasLen;
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

  int32   pageWidth = to - from + 1;

  char   *sequence = new char [length() + 1];
  char   *quality  = new char [length() + 1];
  uint32  len      = 0;

  getConsensus(abacus, sequence, quality, len, length() + 1);

  assert(len == length());


  uint32 numSeqs = abacus->numberOfSequences();

  abSeqBeadIterator     **fit       = new abSeqBeadIterator * [numSeqs];
  uint32                 *fid       = new uint32              [numSeqs];
  char                   *type      = new char                [numSeqs];  //  always 'r' for read
  uint32                 *bgn       = new uint32              [numSeqs];  //  former positions
  uint32                 *end       = new uint32              [numSeqs];



  for (uint32 i=0; i<numSeqs; i++) {
    abSequence  *seq = abacus->getSequence(i);

    assert(seq->multiAlignID() == ident());

    fid[i]  = seq->gkpIdent();
    fit[i]  = NULL;
    type[i] = (seq->isRead() == true) ? 'r' : '?';

    abColumn   *bgnColumn = abacus->getColumn(seq->firstBead());
    abColumn   *endColumn = abacus->getColumn(seq->lastBead());

    bgn[i] = bgnColumn->position();
    bgn[i] = endColumn->position();
  }

  uint32 window_start = from;

  fprintf(F,"\n\n================  MultiAlignment ID %d ==================\n\n", ident().get());

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

            fit[i] = abacus->createSeqBeadIterator(i);
            fit[i]->advance(wi);

            abBeadID bid = fit[i]->next();

            if (bid.isValid()) {
              char pc = abacus->getBase(bid);
              
              if (pc == sequence[wi])
                pc = tolower(pc);
              else
                pc = toupper(pc);

              fprintf(F, "%c", pc);
            }

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
        }


        //  This really looks like a bug.  If we had no iterator, and we started in the middle of
        //  the sequence (HOW?)  we'd print both a base above, and a base below.


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

          } else {
            fprintf(F," ");
            delete fit[i];
            fit[i] = NULL;
          }
        }


        if (wi == window_start + pageWidth - 1)
          fprintf(F," <<< %d (%c)\n", fid[i], type[i]);
      }
    }

    window_start += pageWidth;
  }

  delete [] fit;
  delete [] fid;
  delete [] type;
  delete [] bgn;
  delete [] end;

  delete [] sequence;
  delete [] quality;
}








#if 0

int
GetMANodePositions(int32        mid,
                   MultiAlignT *ma) {

  MANode *manode = GetMANode(manodeStore, mid);

  int32  n_frags   = 0;
  int32  n_unitigs = 0;

  ResetVA_int32(ma->fdelta);
  ResetVA_int32(ma->udelta);

  for (uint32 i=0; i<GetNumFragments(fragmentStore); i++) {
    Fragment *fragment = GetFragment(fragmentStore, i);

    //fprintf(stderr, "GetMANodePositions()--  frag %d ident %d deleted %d\n",
    //        i, fragment->iid, fragment->deleted);

    if (fragment->deleted)
      //  Ejected by contig consensus??
      continue;

    if (fragment->manode == -1)
      //  Fragment not placed in this multialign
      continue;

    assert(fragment->manode == mid);

    int32 bgn = GetColumn(columnStore, (GetBead(beadStore,fragment->firstbead.get()                       ))->column_index)->ma_index;
    int32 end = GetColumn(columnStore, (GetBead(beadStore,fragment->firstbead.get() + fragment->length - 1))->column_index)->ma_index + 1;

    if (fragment->type == AS_READ) {

      //  Not valid for unitig consensus.
      if (fragmentMap) {
        if (FALSE == ExistsInHashTable_AS (fragmentMap, fragment->iid, 0))
          //  Fragment is not in the contig f_list; is in a surrogate.
          continue;

        if (1 != LookupValueInHashTable_AS (fragmentMap, fragment->iid, 0))
          //  Attempting to place a surrogate fragment more than once.
          continue;

        //  Indicate we've placed the fragment.
        ReplaceInHashTable_AS(fragmentMap, fragment->iid, 0, 2, 0);
      }

      IntMultiPos *imp = GetIntMultiPos(ma->f_list, n_frags++);

      imp->ident        = fragment->iid;
      imp->type         = fragment->type;
      imp->position.bgn = (fragment->complement) ? end : bgn;
      imp->position.end = (fragment->complement) ? bgn : end;
      imp->delta        = NULL;
      imp->delta_length = GetFragmentDeltas(i, ma->fdelta, fragment->length);
    }

    if (fragment->type  == AS_UNITIG) {
      IntUnitigPos  *iup = GetIntUnitigPos(ma->u_list, n_unitigs++);

      assert(iup->ident == fragment->iid);

      iup->position.bgn = (fragment->complement) ? end : bgn;
      iup->position.end = (fragment->complement) ? bgn : end;
      iup->delta        = NULL;
      iup->delta_length = GetFragmentDeltas(i, ma->udelta, fragment->length);
    }
  }

  //  Because contig consensus might have ejected fragments that don't align, the new list can be
  //  shorter than the original list.
  //
  ResetToRangeVA_IntMultiPos(ma->f_list, n_frags);

  //  Set delta pointers into the VA.

  int32  fdeltapos = 0;
  int32  udeltapos = 0;

  n_frags   = 0;
  n_unitigs = 0;

  for (uint32 i=0; i<GetNumFragments(fragmentStore); i++) {
    Fragment *fragment = GetFragment(fragmentStore, i);

    if (fragment->deleted)
      //  Ejected by contig consensus??
      continue;

    if (fragment->manode == -1)
      //  Fragment not placed in this multialign
      continue;

    if (fragment->type == AS_READ) {

      //  Not valid for unitig consensus

      if (fragmentMap) {
        if (FALSE == ExistsInHashTable_AS(fragmentMap, fragment->iid, 0))
          continue;

        // all of the contig's fragments should've had their value set to 2 in previous block

        assert(2 == LookupValueInHashTable_AS(fragmentMap, fragment->iid, 0));
        DeleteFromHashTable_AS(fragmentMap, fragment->iid, 0);
      }

      IntMultiPos *imp = GetIntMultiPos(ma->f_list, n_frags++);

      imp->delta = (imp->delta_length == 0) ? NULL : Getint32(ma->fdelta, fdeltapos);

      fdeltapos += imp->delta_length;
    }

    if (fragment->type == AS_UNITIG) {
      IntUnitigPos  *iup = GetIntUnitigPos(ma->u_list, n_unitigs++);

      iup->delta = (iup->delta_length == 0) ? NULL : Getint32(ma->udelta, udeltapos);

      udeltapos += iup->delta_length;
    }
  }
}



#endif
