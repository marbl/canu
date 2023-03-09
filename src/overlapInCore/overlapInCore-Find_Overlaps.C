
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

//  Add information for the match in  ref  to the list
//  starting at subscript  (* start). The matching window begins
//  offset  bytes from the beginning of this string.


static
void
workArea::addMatch(kmerRef ref,
                   int * start,
                   int offset,             //  Should be increasing; position in read we're processing
                   int * consistent,
                   Work_Area_t * WA) {

  int  * p, save;
  int  diag = 0, new_diag, expected_start = 0, num_checked = 0;
  int  move_to_front = false;

  new_diag = ref._stringPos - offset;  //  diagonal in hash_read vs query_read alignmment

  // walk through all the entries for this read
  // 
  for (p = start;  (* p) != 0;  p = & (WA->Match_Node_Space [(* p)].Next)) {

    //  compute the expected start of an extension of this match
    //  and the diagonal we're on
    expected_start = WA->Match_Node_Space [(* p)].Start  + WA->Match_Node_Space [(* p)].Len - G.Kmer_Len + 1;
    diag           = WA->Match_Node_Space [(* p)].Offset - WA->Match_Node_Space [(* p)].Start;

    //  if the expected start is before the start of the match, stop.
    //  why??
    if (expected_start < offset)
      break;

    //  if the expected start is exactly the start of the new match
    if (expected_start == offset) {

      //  and it's on the same diagonal,
      //    increase the length of this match by one
      //    possibly move this entry to the start of the list
      //    and return
      if (new_diag == diag) {
        WA->Match_Node_Space [(* p)].Len += 1;

        //  if we ...
        //    save the index we're looking at
        //    move to the next one
        //    set the Next pointer to point to the 'start' of the list
        //    and make the start of the list be the index we're looking at
        //  so indeed we just move this Node_Space to the start of the list
        //
        if (move_to_front) {
          save = (* p);
          (* p) = WA->Match_Node_Space [(* p)].Next;
          WA->Match_Node_Space [save].Next = (* start);
          (* start) = save;
        }
        return;
      }

      //  but if on a different diagonal
      //    flag it so we move it to the start of the list if we ever extend it
      //    i suspect this keeps the Match_Node_Space list sorted by decreasing expected_start
      else
        move_to_front = true;
    }

    num_checked ++;
  }

  //  if here we either failed to find an exact 'expected_start' match on the correct diagonal
  //  or all the rest of the Match_Node_Space entries end before we start

  // Allocate more match node space if needed
  if (WA->Next_Avail_Match_Node == WA->Match_Node_Size) {
    int32          newSize  = WA->Match_Node_Size * 2;
    Match_Node_t  *newSpace = new Match_Node_t [newSize];

    memcpy(newSpace, WA->Match_Node_Space, sizeof(Match_Node_t) * WA->Match_Node_Size);

    delete [] WA->Match_Node_Space;

    WA->Match_Node_Size  = newSize;
    WA->Match_Node_Space = newSpace;
  }

  //  if a bunch of stuff is true set non-consistent
  //    start isn't the first node in the array
  //      we checked at least one node -- the first node started after this kmer
  //      a very different diagonal
  //      a gap in the kmers
  //  then not consistent

  if ((* start) != 0
      && (num_checked > 0
          || abs (diag - new_diag) > 3
          || offset < expected_start + G.Kmer_Len - 2))
    (* consistent) = false;

  // remember the address of the current head, set 'start' to the next free node
  save = (* start);
  (* start) = WA->Next_Avail_Match_Node;
  WA->Next_Avail_Match_Node ++;

  // populate the free node
  WA->Match_Node_Space [(* start)].Offset = ref._stringPos;  //  start of align in hash table read
  WA->Match_Node_Space [(* start)].Len    = G.Kmer_Len;      //  length of hit
  WA->Match_Node_Space [(* start)].Start  = offset;          //  start of align in query read
  WA->Match_Node_Space [(* start)].Next   = save;            //  pointer to the rest of the hits

#if 0
  fprintf(stderr, "addMatch()-- %3d offset %d len %d start %d next %d\n",
          *start,
          WA->Match_Node_Space [(* start)].Offset,
          WA->Match_Node_Space [(* start)].Len,
          WA->Match_Node_Space [(* start)].Start,
          WA->Match_Node_Space [(* start)].Next);
#endif
}



//  Add information for Ref and all its matches to the global hash table in String_Olap_Space. Grow
//  the space if necessary. The matching window begins Offset bytes from the beginning of this
//  string.

void
workArea::addReference(kmerRef kmer,
                       uint64  offset) {

  uint32  Prev;

  uint32 StrNum = kmer._stringID;
  uint32 Sub    = (StrNum ^ (StrNum >> STRING_OLAP_SHIFT)) & STRING_OLAP_MASK;

  while ((WA->String_Olap_Space[Sub].Full) &&
         (WA->String_Olap_Space[Sub].String_Num != StrNum)) {

    Prev = Sub;
    Sub = WA->String_Olap_Space[Sub].Next;

    if (Sub == 0) {
      if (WA->Next_Avail_String_Olap == WA->String_Olap_Size) {
        int32          newSize  = WA->String_Olap_Size * 2;
        String_Olap_t *newSpace = new String_Olap_t [newSize];

        memcpy(newSpace, WA->String_Olap_Space, sizeof(String_Olap_t) * WA->String_Olap_Size);

        delete [] WA->String_Olap_Space;

        WA->String_Olap_Size  = newSize;
        WA->String_Olap_Space = newSpace;
      }

      Sub = WA->Next_Avail_String_Olap++;

      WA->String_Olap_Space [Prev].Next = Sub;
      WA->String_Olap_Space[Sub].Full = false;

      break;
    }
  }

  if (! WA->String_Olap_Space[Sub].Full) {
    WA->String_Olap_Space[Sub].String_Num = StrNum;
    WA->String_Olap_Space[Sub].Match_List = 0;
    WA->String_Olap_Space[Sub].diag_sum = 0.0;
    WA->String_Olap_Space[Sub].diag_ct = 0;
    WA->String_Olap_Space[Sub].diag_bgn = AS_MAX_READLEN;
    WA->String_Olap_Space[Sub].diag_end = 0;
    WA->String_Olap_Space[Sub].Next = 0;
    WA->String_Olap_Space[Sub].Full = true;
    WA->String_Olap_Space[Sub].consistent = true;
  }

  int consistent = WA->String_Olap_Space[Sub].consistent;

  WA->String_Olap_Space[Sub].diag_sum += (double)kmer._stringPos - offset;
  WA->String_Olap_Space[Sub].diag_ct++;

  if (WA->String_Olap_Space[Sub].diag_bgn > offset) WA->String_Olap_Space[Sub].diag_bgn = offset;
  if (WA->String_Olap_Space[Sub].diag_end < offset) WA->String_Olap_Space[Sub].diag_end = offset;

  //                start of list of exact matches in match node space
  //
  addMatch(kmer, & (WA->String_Olap_Space[Sub].Match_List), offset, & consistent, WA);

  WA->String_Olap_Space[Sub].consistent = consistent;

  return;
}




//  Search for string  S  with hash key  Key  in the global
//  Hash_Table  starting at subscript  Sub. Return the matching
//  reference in the hash table if there is one, or else a reference
//  with the  Empty bit set true.  Set  (* Where)  to the subscript in
//  Extra_Ref_Space  where the reference was found if it was found there.
//  Set  (* hi_hits)  to  true  if hash table entry is found but is empty
//  because it was screened out, otherwise set to false.

kmerRef
Hash_Find(uint64 key, char *s) {
  uint64  bucket   = Hash_Bucket_t::hashFunction(key);
  uint64  probe    = Hash_Bucket_t::hashProbe(key);
  uint8   keycheck = Hash_Bucket_t::keyCheck(key);

  //  If the check bit isn't set for this bucket, then the kmer isn't in this
  //  bucket or any of it's probe buckets.

  if (Hash_Table[bucket].isNotPresent(key) == true)
    return(kmerRef());

  //  Iterate until we find the kmer, or we've visited every bucket in the
  //  table and decide the table is full.

  for (int64 Ct=0; Ct < HASH_TABLE_SIZE; Ct++) {
    for (uint32 ee=0; ee < Hash_Table[bucket]._bucketLen; ee++) {
      kmerRef  hashref = Hash_Table[bucket]._bucket[ee];
      char         *t       = basesData + hashref.position();

      if ((Hash_Table[bucket]_bcheck[ee] != keycheck) ||   //  If check fails, kmer can't be here,
          (strncmp(s, t, G.Kmer_Len) != 0))               //  otherwise verify it's the same kmer.
        continue;

      return(hashref);
    }

    //  If the bucket isn't full, we're done, no need to probe ahead to the
    //  next bucket.  Otherwise, probe ahead!

    if (Hash_Table[bucket]._bucketLen == ENTRIES_PER_BUCKET)
      return(kmerRef());

    bucket = (bucket + probe) % HASH_TABLE_SIZE;
  }

  //  We checked everywhere we could and didn't find it.  Bummer, eh.

  return(kmerRef());
}




//  Lookup each kmer in the hash table.
//
//  If the kmer is frequent and we're near the end of the read, flag the
//  read as such and ignore the hits.
//
//  If it isn't frequent, zip down the linear reference array adding hits
//  for each position.
//
void
workArea::findOverlaps(char       *readSeq,
                       uint32      readLen,
                       uint32      readID,
                       bool        doForward,
                       hashTable  *HT) {

  memset(String_Olap_Space, 0, STRING_OLAP_MODULUS * sizeof (String_Olap_t));

  Next_Avail_String_Olap = STRING_OLAP_MODULUS;
  Next_Avail_Match_Node  = 1;

  left_end_screened  = false;
  right_end_screened = false;

  A_Olaps_For_Frag = 0;
  B_Olaps_For_Frag = 0;

  keyIterator    it(readSeq, HT->kmerSize());

  while (it.next() == true) {
    kmerRef  hashref = Hash_Find(it.key, it.merBgn);
    uint64   refiter = hashref.getLinearPointer();

    if (hashref._ignore) {
      int32  bgnSize =            it.ref._stringPos;
      int32  endSize = Frag_Len - it.ref._stringPos - G.Kmer_Len + 1;

      if (bgnSize < HOPELESS_MATCH)   left_end_screened  = true;
      if (endSize < HOPELESS_MATCH)   right_end_screened = true;

      continue;
    }

    hashref = Extra_Ref_Space[refiter++];

    while (hashref._isLast == false) {
      if (Frag_Num < HT->_bgnRead + hashref._stringID)
        addReference(hashref,                            //  the hit in the hash table
                     it.ref._stringPos - G.Kmer_Len);    //  start position of kmer in the read

      hashref = Extra_Ref_Space[refiter++];
    }
  }
}
