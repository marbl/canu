#include "abAbacus.H"


#if 0

//  If we're adding a unitig (with reads in it) drill into the unitig and add each
//  of its reads (at least, I think that's what this is doing)

static
int
abAbacus::SetUngappedFragmentPositions(FragType type,int32 n_frags, MultiAlignT *uma) {

  int32 num_frags   = GetNumIntMultiPoss(uma->f_list);
  int32 num_unitigs = GetNumIntUnitigPoss(uma->u_list);

  HashTable_AS *unitigFrags = CreateScalarHashTable_AS();

  int32 num_columns   = GetMultiAlignLength(uma);
  int32 ungapped_pos  = 0;

  int32 *gapped_positions = new int32 [num_columns + 1];
  char  *consensus        = Getchar(uma->consensus,0);

  for (int32 i=0; i<num_columns+1; i++) {
    gapped_positions[i] = ungapped_pos;

    if (consensus[i] != '-')
      ungapped_pos++;
  }

  //  Remember the first fragment we add.
  int32 first_frag = GetNumCNS_AlignedContigElements(fragment_positions);

  for (int32 ifrag=0; ifrag<num_frags; ifrag++) {
    CNS_AlignedContigElement epos;
    IntMultiPos *frag = GetIntMultiPos(uma->f_list, ifrag);

    if (ExistsInHashTable_AS(unitigFrags, frag->ident, 0)) {
      fprintf(stderr,"SetUngappedFragmentPositions()-- ident %d already in hashtable\n", frag->ident);
      assert(0);
    }
    if (HASH_SUCCESS != InsertInHashTable_AS(unitigFrags, frag->ident, 0, 1, 0)) {
      fprintf(stderr,"SetUngappedFragmentPositions()-- Failure to insert ident %d in hashtable\n", frag->ident);
      assert(0);
    }

    assert(frag->position.bgn >= 0);
    assert(frag->position.bgn < num_columns + 1);
    assert(frag->position.end >= 0);
    assert(frag->position.end < num_columns + 1);

    epos.frg_or_utg                  = CNS_ELEMENT_IS_FRAGMENT;
    epos.idx.fragment.frgIdent       = frag->ident;
    epos.idx.fragment.frgType        = frag->type;
    epos.idx.fragment.frgContained   = frag->contained;
    epos.idx.fragment.frgInUnitig    = (type == AS_CONTIG) ? -1 : uma->maID;
    epos.position.bgn                = gapped_positions[frag->position.bgn];
    epos.position.end                = gapped_positions[frag->position.end];

    //fprintf(stderr, "SetUngappedFragmentPositions()-- FRG id=%d type=%c pos=%d,%d (orig pos=%d,%d)\n",
    //        frag->ident, frag->type, epos.position.bgn, epos.position.end, frag->position.bgn, frag->position.end);

    //  Adjust the ungapped position if we fall within a gap
    //
    if (epos.position.bgn == epos.position.end) {
      fprintf(stderr,"SetUngappedFragmentPositions()-- Encountered bgn==end=="F_S32" in ungapped coords within SetUngappedFragmentPositions for "F_CID "(gapped coords "F_S32","F_S32")\n",
              epos.position.bgn,frag->ident,frag->position.bgn,frag->position.end);
      assert(frag->position.bgn != frag->position.end);

      if (frag->position.bgn < frag->position.end) {
        if (epos.position.bgn > 0)
          epos.position.bgn--;
        else
          epos.position.end++;
      } else {
        if (epos.position.end > 0)
          epos.position.end--;
        else
          epos.position.bgn++;
      }
      fprintf(stderr,"SetUngappedFragmentPositions()--   Reset to "F_S32","F_S32"\n",
              epos.position.bgn,
              epos.position.end);
    }

    AppendVA_CNS_AlignedContigElement(fragment_positions, &epos);
  }


  for (int32 ifrag=0; ifrag < num_unitigs; ifrag++){
    CNS_AlignedContigElement epos;
    IntUnitigPos *unitig = GetIntUnitigPos(uma->u_list, ifrag);

    epos.frg_or_utg           = CNS_ELEMENT_IS_UNITIG;
    epos.idx.unitig.utgIdent  = unitig->ident;
    epos.idx.unitig.utgType   = unitig->type;
    epos.position.bgn         = gapped_positions[unitig->position.bgn];
    epos.position.end         = gapped_positions[unitig->position.end];

    //fprintf(stderr, "SetUngappedFragmentPositions()-- UTG id=%d type=%c pos=%d,%d (orig pos=%d,%d)\n",
    //        unitig->ident, unitig->type, epos.position.bgn, epos.position.end, unitig->position.bgn, unitig->position.end);

    AppendVA_CNS_AlignedContigElement(fragment_positions,&epos);
  }

  //  This is used only by ReplaceEndUnitigInContig().  Mark fragments in the "anchoring" contig
  //  that belong to this unitig.
  //
  if (type != AS_CONTIG) {
    Fragment *anchor = GetFragment(fragmentStore,0);

    if ((anchor != NULL) &&
        (anchor->type == AS_CONTIG)) {
      CNS_AlignedContigElement *af = GetCNS_AlignedContigElement(fragment_positions, anchor->components);

      for (int32 ifrag=0; ifrag < anchor->n_components; ifrag++, af++) {
        if ((af->frg_or_utg == CNS_ELEMENT_IS_FRAGMENT) &&
            (ExistsInHashTable_AS(unitigFrags, af->idx.fragment.frgIdent, 0)))
          af->idx.fragment.frgInUnitig = uma->maID;
      }
    }
  }

  DeleteHashTable_AS(unitigFrags);
  delete [] gapped_positions;

  return first_frag;
}



int32
abAbacus::addContig() {

      if (tigStore)
        uma = tigStore->loadMultiAlign(iid, type == AS_UNITIG);
      if (uma == NULL)
        fprintf(stderr,"Lookup failure in CNS: MultiAlign for unitig %d could not be found.\n",iid);
      assert(uma != NULL);

      //  Contigs used to be added gapped, unitigs as ungapped.
      //  This caused no end of trouble in MergeMultiAligns and
      //  ReplaceEndUnitigInContig.

      ResetVA_char(ungappedSequence);
      ResetVA_char(ungappedQuality);

      GetMultiAlignUngappedConsensus(uma, ungappedSequence, ungappedQuality);

      sequence = Getchar(ungappedSequence,0);
      quality = Getchar(ungappedQuality,0);

      fragment.length = GetMultiAlignUngappedLength(uma);

      fragment.utype = (type == AS_UNITIG) ? utype : AS_OTHER_UNITIG;

      fragment.n_components = GetNumIntMultiPoss(uma->f_list) + GetNumIntUnitigPoss(uma->u_list);
      fragment.components   = SetUngappedFragmentPositions(type, fragment.n_components, uma);

      //fprintf(stderr, "AppendFragToLocalStore()-- TIG %d len=%d\n", iid, fragment.length);
}

#endif





abSeqID
abAbacus::addRead(gkStore *gkpStore,
                  uint32   readID,
                  uint32   askip, uint32 bskip,
                  bool     complemented) {

  //  Grab the read

  gkRead      *read     = gkpStore->gkStore_getRead(readID);
  gkReadData  *readData = new gkReadData;

  //fprintf(stderr, "abAbacus::addRead()--  want readID=%u, store returned read %u\n",
  //        readID, read->gkRead_readID());

  gkpStore->gkStore_loadReadData(read, readData);

  uint32  seqLen = read->gkRead_sequenceLength() - askip - bskip;

  //  Tell abacus about it.

  increaseArray(_sequences, _sequencesLen, _sequencesMax, 1);

  abSequence   *ns = _sequences + _sequencesLen;

  ns->initialize(readID, _sequencesLen++, seqLen, complemented, _basesLen, _beadsLen);

  //  Make a complement table

  char inv[256] = {0};

  inv['a'] = 't';  inv['A'] = 'T';
  inv['c'] = 'g';  inv['C'] = 'G';
  inv['g'] = 'c';  inv['G'] = 'C';
  inv['t'] = 'a';  inv['T'] = 'A';
  inv['n'] = 'n';  inv['N'] = 'N';
  inv['-'] = '-';

  //  Stash the bases/quals

  {
    char  *seq = readData->gkReadData_getSequence()  + ((complemented == false) ? askip : bskip);
    char  *qlt = readData->gkReadData_getQualities() + ((complemented == false) ? askip : bskip);

    while (_basesMax <= _basesLen + seqLen + 1)
      resizeArrayPair(_bases, _quals, _basesLen, _basesMax, 2 * _basesMax);

    for (uint32 ii=0; ii<seqLen; ii++)
      assert((seq[ii] == 'A') ||
             (seq[ii] == 'C') ||
             (seq[ii] == 'G') ||
             (seq[ii] == 'T') ||
             (seq[ii] == 'N'));

    if (complemented == false)
      for (uint32 ii=0, pp=_basesLen; ii<seqLen; ii++, pp++, _basesLen++) {
        _bases[pp] = seq[ii];
        _quals[pp] = qlt[ii];

        assert(CNS_MIN_QV + '0' <= _quals[pp]);
        assert(_quals[pp] <= CNS_MAX_QV + '0');
      }

    else
      for (uint32 ii=seqLen, pp=_basesLen; ii-->0; pp++, _basesLen++) {
        _bases[pp] = inv[ seq[ii] ];
        _quals[pp] =      qlt[ii];

        assert(CNS_MIN_QV + '0' <= _quals[pp]);
        assert(_quals[pp] <= CNS_MAX_QV + '0');
      }

    _bases[_basesLen] = 0;  //  NUL terminate the strings so we can use them in aligners
    _quals[_basesLen] = 0;
    _basesLen++;
  }

  delete readData;

  //  Make beads for each base, set the pointer to the first bead

  {
    increaseArray(_beads, _beadsLen, _beadsMax, seqLen);

    uint32  firstBead = _beadsLen;

    for (uint32 bp=0; bp<ns->length(); bp++, _beadsLen++) {
      _beads[_beadsLen].boffset.set(_beadsLen);                   //  Offset into the beads array
      _beads[_beadsLen].soffset.set(ns->firstBase().get() + bp);  //  Offset into the sequence array
      _beads[_beadsLen].foffset = bp;                             //  Offset into the read itself

      //  Check that nothing odd happened with the ident.

      assert(_beads[_beadsLen].ident().get() == ns->firstBead().get() + bp);

      //  Set previous/next bead appropriately.

      if (_beads[_beadsLen].foffset == 0)
        _beads[_beadsLen].prev = abBeadID();
      else
        _beads[_beadsLen].prev.set(_beads[_beadsLen].ident().get() - 1);

      if (_beads[_beadsLen].foffset == ns->length() - 1)
        _beads[_beadsLen].next = abBeadID();
      else
        _beads[_beadsLen].next.set(_beads[_beadsLen].ident().get() + 1);

      _beads[_beadsLen].up           = abBeadID();   //  No up bead yet.
      _beads[_beadsLen].down         = abBeadID();   //  No down bead yet.

      _beads[_beadsLen].frag_index   = ns->ident();  //  Bead is for this read idx.
      _beads[_beadsLen].column_index = abColID();    //  Isn't in a column yet.
    }
  }

  assert(_beads[_beadsLen-1].ident() == ns->lastBead());

  //  Return the (internal) index we saved this read at.

#if 0
  fprintf(stderr, "read %d firstBead %d lastBead %d _basesLen %u\n",
          ns->ident().get(),
          ns->firstBead().get(),
          ns->lastBead().get(),
          _basesLen);
#endif

  return(ns->ident());
}


