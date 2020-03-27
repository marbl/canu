
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

#include "overlapInCore.H"
#include "sequence.H"

//  Find and output all overlaps between strings in store and those in the global hash table.
//  This is the entry point for each compute thread.

void *
Process_Overlaps(void *ptr){
  Work_Area_t  *WA = (Work_Area_t *)ptr;

  uint32        seqptrLen = 0;
  uint32        seqptrMax = AS_MAX_READLEN + 1;
  char         *seqptr    = new char [seqptrMax];
  char         *bases     = new char [AS_MAX_READLEN + 1];

  while (WA->bgnID < G.endRefID) {
    WA->overlapsLen                = 0;

    WA->Total_Overlaps             = 0;
    WA->Contained_Overlap_Ct       = 0;
    WA->Dovetail_Overlap_Ct        = 0;

    WA->Kmer_Hits_Without_Olap_Ct  = 0;
    WA->Kmer_Hits_With_Olap_Ct     = 0;
    WA->Kmer_Hits_Skipped_Ct       = 0;
    WA->Multi_Overlap_Ct           = 0;

    fprintf(stderr, "Thread %02u processes reads " F_U32 "-" F_U32 "\n",
            WA->thread_id, WA->bgnID, WA->endID);

    for (uint32 fi=WA->bgnID; fi<=WA->endID; fi++) {
      uint32  libID   = WA->readStore->sqStore_getLibraryIDForRead(fi);
      uint32  readLen = WA->readCache->sqCache_getLength(fi);

      //  Load sequence/quality data
      //  Duplicated in Build_Hash_Index()

      if ((libID < G.minLibToRef) ||
          (libID > G.maxLibToRef))
        continue;

      if (readLen < G.Min_Olap_Len)
        continue;

      WA->readCache->sqCache_getSequence(fi, seqptr, seqptrLen, seqptrMax);

      for (uint32 i=0; i<readLen; i++)
        bases[i] = tolower(seqptr[i]);

      bases[readLen] = 0;

      assert(strlen(bases) == readLen);

      //  Generate overlaps.

      Find_Overlaps(bases, readLen, fi, FORWARD, WA);

      reverseComplementSequence(bases, readLen);

      Find_Overlaps(bases, readLen, fi, REVERSE, WA);
    }

    //  Write out this block of overlaps, no need to keep them in core!
    //  While we have a mutex, also find the next block of things to process.

    fprintf(stderr, "Thread %02u writes    reads " F_U32 "-" F_U32 " (" F_U64 " overlaps " F_U64 "/" F_U64 "/" F_U64 " kmer hits with/without overlap/skipped)\n",
            WA->thread_id, WA->bgnID, WA->endID,
            WA->overlapsLen,
            WA->Kmer_Hits_With_Olap_Ct, WA->Kmer_Hits_Without_Olap_Ct, WA->Kmer_Hits_Skipped_Ct);

    //  Flush any remaining overlaps and update statistics.

#pragma omp critical
    {
      for (int zz=0; zz<WA->overlapsLen; zz++)
        Out_BOF->writeOverlap(WA->overlaps + zz);

      WA->overlapsLen = 0;

      Total_Overlaps            += WA->Total_Overlaps;
      Contained_Overlap_Ct      += WA->Contained_Overlap_Ct;
      Dovetail_Overlap_Ct       += WA->Dovetail_Overlap_Ct;

      Kmer_Hits_Without_Olap_Ct += WA->Kmer_Hits_Without_Olap_Ct;
      Kmer_Hits_With_Olap_Ct    += WA->Kmer_Hits_With_Olap_Ct;
      Kmer_Hits_Skipped_Ct      += WA->Kmer_Hits_Skipped_Ct;
      Multi_Overlap_Ct          += WA->Multi_Overlap_Ct;

      WA->bgnID = G.curRefID;
      WA->endID = G.curRefID + G.perThread - 1;

      if (WA->endID > G.endRefID)
        WA->endID = G.endRefID;

      G.curRefID = WA->endID + 1;
    }
  }

  delete [] bases;
  delete [] seqptr;

  return(ptr);
}


