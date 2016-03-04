
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_OVL/AS_OVL_overlap_common.h
 *    src/AS_OVM/overlapInCore-Build_Hash_Index.C
 *
 *  Modifications by:
 *
 *    Michael Schatz on 2004-SEP-23
 *      are Copyright 2004 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2005-MAR-22
 *      are Copyright 2005 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUN-16 to 2013-AUG-01
 *      are Copyright 2005-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2005-JUL-15 to 2007-NOV-20
 *      are Copyright 2005,2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2006-MAR-27 to 2006-AUG-21
 *      are Copyright 2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Art Delcher on 2007-FEB-13
 *      are Copyright 2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2007-AUG-27 to 2009-JAN-16
 *      are Copyright 2007,2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2011-MAR-08 to 2011-JUN-17
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-DEC-15 to 2015-AUG-25
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "overlapInCore.H"

#include "AS_UTL_reverseComplement.H"



//  Add string  s  as an extra hash table string and return
//  a single reference to the beginning of it.
static
String_Ref_t
Add_Extra_Hash_String(const char *s) {
  String_Ref_t  ref = 0;
  String_Ref_t  sub = 0;

  int  len;

  uint32 new_len = Used_Data_Len + G.Kmer_Len;

  if (Extra_String_Subcount < MAX_EXTRA_SUBCOUNT) {
    sub = String_Ct + Extra_String_Ct - 1;

  } else {
    sub = String_Ct + Extra_String_Ct;

    if (sub >= String_Start_Size) {
      uint64  n = max(sub * 1.1, String_Start_Size * 1.5);

      //fprintf(stderr, "REALLOC String_Start from "F_U64" to "F_U64"\n", String_Start_Size, n);
      resizeArray(String_Start, String_Start_Size, String_Start_Size, n);
    }

    String_Start[sub] = Used_Data_Len;

    Extra_String_Ct++;
    Extra_String_Subcount = 0;
    new_len++;
  }

  if (new_len >= Extra_Data_Len) {
    uint64  n = max(new_len * 1.1, Extra_Data_Len * 1.5);

    //fprintf(stderr, "REALLOC basesData from "F_U64" to "F_U64"\n", Extra_Data_Len, n);
    resizeArray(basesData, Extra_Data_Len, Extra_Data_Len, n);
  }

  strncpy(basesData + String_Start[sub] + G.Kmer_Len * Extra_String_Subcount, s, G.Kmer_Len + 1);

  Used_Data_Len = new_len;

  setStringRefStringNum(ref, sub);

  if (sub > MAX_STRING_NUM) {
    fprintf(stderr, "Too many skip kmer strings for hash table.\n");
    fprintf(stderr, "Try skipping hopeless check (-z option)\n");
    fprintf(stderr, "Exiting\n");
    exit (1);
  }

  setStringRefOffset(ref, (String_Ref_t)Extra_String_Subcount * (String_Ref_t)G.Kmer_Len);

  assert(Extra_String_Subcount * G.Kmer_Len < OFFSET_MASK);

  setStringRefLast(ref,  (uint64)1);
  setStringRefEmpty(ref, TRUELY_ONE);

  Extra_String_Subcount++;

  return(ref);
}





//  Mark  left/right_end_screened in global  String_Info  for
//   ref  and everything in its list, if they occur near
//  enough to the end of the string.

static
void
Mark_Screened_Ends_Single(String_Ref_t ref) {
  int32 s_num = getStringRefStringNum(ref);
  int32 len = String_Info[s_num].length;

  if (getStringRefOffset(ref) < HOPELESS_MATCH)
    String_Info[s_num].lfrag_end_screened = TRUE;

  if (len - getStringRefOffset(ref) - G.Kmer_Len + 1 < HOPELESS_MATCH)
    String_Info[s_num].rfrag_end_screened = TRUE;
}



static
void
Mark_Screened_Ends_Chain(String_Ref_t ref) {

  Mark_Screened_Ends_Single (ref);

  while (! getStringRefLast(ref)) {
    ref = nextRef[(String_Start[getStringRefStringNum(ref)] + getStringRefOffset(ref)) / (HASH_KMER_SKIP + 1)];
    Mark_Screened_Ends_Single (ref);
  }
}


//  Set the  empty  bit to true for the hash table entry
//  corresponding to string  s  whose hash key is  key .
//  Also set global  String_Info.left/right_end_screened
//  true if the entry occurs near the left/right end, resp.,
//  of the string in the hash table.  If not found, add an
//  entry to the hash table and mark it empty.
static
void
Hash_Mark_Empty(uint64 key, char * s) {
  String_Ref_t  h_ref;
  char  * t;
  unsigned char  key_check;
  int64  ct, probe;
  int64  sub;
  int  i, shift;

  sub = HASH_FUNCTION (key);
  key_check = KEY_CHECK_FUNCTION (key);
  probe = PROBE_FUNCTION (key);

  ct = 0;
  do {
    for (i = 0;  i < Hash_Table[sub].Entry_Ct;  i ++)
      if (Hash_Table[sub].Check[i] == key_check) {
        h_ref = Hash_Table[sub].Entry[i];
        t = basesData + String_Start[getStringRefStringNum(h_ref)] + getStringRefOffset(h_ref);
        if (strncmp (s, t, G.Kmer_Len) == 0) {
          if (! getStringRefEmpty(Hash_Table[sub].Entry[i]))
            Mark_Screened_Ends_Chain (Hash_Table[sub].Entry[i]);
          setStringRefEmpty(Hash_Table[sub].Entry[i], TRUELY_ONE);
          return;
        }
      }
    assert (i == Hash_Table[sub].Entry_Ct);
    if (Hash_Table[sub].Entry_Ct < ENTRIES_PER_BUCKET) {
      // Not found
      if (G.Use_Hopeless_Check) {
        Hash_Table[sub].Entry[i] = Add_Extra_Hash_String (s);
        setStringRefEmpty(Hash_Table[sub].Entry[i], TRUELY_ONE);
        Hash_Table[sub].Check[i] = key_check;
        Hash_Table[sub].Entry_Ct ++;
        Hash_Table[sub].Hits[i] = 0;
        Hash_Entries ++;
        shift = HASH_CHECK_FUNCTION (key);
        Hash_Check_Array[sub] |= (((Check_Vector_t) 1) << shift);
      }
      return;
    }
    sub = (sub + probe) % HASH_TABLE_SIZE;
  }  while (++ ct < HASH_TABLE_SIZE);

  fprintf (stderr, "ERROR:  Hash table full\n");
  assert (FALSE);
}



//  Set  Empty  bit true for all entries in global  Hash_Table
//  that match a kmer in file  Kmer_Skip_File .
//  Add the entry (and then mark it empty) if it's not in  Hash_Table.
static
void
Mark_Skip_Kmers(void) {
  uint64  key;
  char  line[MAX_LINE_LEN];
  int  ct = 0;

  rewind (G.Kmer_Skip_File);

  while (fgets (line, MAX_LINE_LEN, G.Kmer_Skip_File) != NULL) {
    int  i, len;

    ct ++;
    len = strlen (line) - 1;
    if (line[0] != '>' || line[len] != '\n') {
      fprintf (stderr, "ERROR:  Bad line %d in kmer skip file\n", ct);
      fputs (line, stderr);
      exit (1);
    }

    if (fgets (line, MAX_LINE_LEN, G.Kmer_Skip_File) == NULL) {
      fprintf (stderr, "ERROR:  Bad line after %d in kmer skip file\n", ct);
      exit (1);
    }
    ct ++;
    len = strlen (line) - 1;
    if (len != G.Kmer_Len || line[len] != '\n') {
      fprintf (stderr, "ERROR:  Bad line %d in kmer skip file\n", ct);
      fputs (line, stderr);
      exit (1);
    }
    line[len] = '\0';

    //if ((ct % 200000) == 0)
    //  fprintf(stderr, "Loaded skip %10d '%s'\n", ct/2, line);

    key = 0;
    for (i = 0;  i < len;  i ++) {
      line[i] = tolower (line[i]);
      key |= (uint64) (Bit_Equivalent[(int) line[i]]) << (2 * i);
    }
    Hash_Mark_Empty (key, line);

    reverseComplementSequence (line, len);
    key = 0;
    for (i = 0;  i < len;  i ++)
      key |= (uint64) (Bit_Equivalent[(int) line[i]]) << (2 * i);
    Hash_Mark_Empty (key, line);
  }

  fprintf (stderr, "String_Ct = "F_U64"  Extra_String_Ct = "F_U64"  Extra_String_Subcount = "F_U64"\n",
           String_Ct, Extra_String_Ct, Extra_String_Subcount);
  fprintf (stderr, "Read %d kmers to mark to skip\n", ct / 2);
}





//  Insert  Ref  with hash key  Key  into global  Hash_Table .
//  Ref  represents string  S .
static
void
Hash_Insert(String_Ref_t Ref, uint64 Key, char * S) {
  String_Ref_t  H_Ref;
  char  * T;
  int  Shift;
  unsigned char  Key_Check;
  int64  Ct, Probe, Sub;
  int  i;

  Sub = HASH_FUNCTION (Key);
  Shift = HASH_CHECK_FUNCTION (Key);
  Hash_Check_Array[Sub] |= (((Check_Vector_t) 1) << Shift);
  Key_Check = KEY_CHECK_FUNCTION (Key);
  Probe = PROBE_FUNCTION (Key);

  Ct = 0;
  do {
    for (i = 0;  i < Hash_Table[Sub].Entry_Ct;  i ++)
      if (Hash_Table[Sub].Check[i] == Key_Check) {
        H_Ref = Hash_Table[Sub].Entry[i];
        T = basesData + String_Start[getStringRefStringNum(H_Ref)] + getStringRefOffset(H_Ref);
        if (strncmp (S, T, G.Kmer_Len) == 0) {
          if (getStringRefLast(H_Ref)) {
            Extra_Ref_Ct ++;
          }
          nextRef[(String_Start[getStringRefStringNum(Ref)] + getStringRefOffset(Ref)) / (HASH_KMER_SKIP + 1)] = H_Ref;
          Extra_Ref_Ct ++;
          setStringRefLast(Ref, TRUELY_ZERO);
          Hash_Table[Sub].Entry[i] = Ref;

          if (Hash_Table[Sub].Hits[i] < HIGHEST_KMER_LIMIT)
            Hash_Table[Sub].Hits[i] ++;

          return;
        }
      }
    if (i != Hash_Table[Sub].Entry_Ct) {
      fprintf (stderr, "i = %d  Sub = " F_S64 "  Entry_Ct = %d\n",
               i, Sub, Hash_Table[Sub].Entry_Ct);
    }
    assert (i == Hash_Table[Sub].Entry_Ct);
    if (Hash_Table[Sub].Entry_Ct < ENTRIES_PER_BUCKET) {
      setStringRefLast(Ref, TRUELY_ONE);
      Hash_Table[Sub].Entry[i] = Ref;
      Hash_Table[Sub].Check[i] = Key_Check;
      Hash_Table[Sub].Entry_Ct ++;
      Hash_Entries ++;
      Hash_Table[Sub].Hits[i] = 1;
      return;
    }
    Sub = (Sub + Probe) % HASH_TABLE_SIZE;
  }  while (++ Ct < HASH_TABLE_SIZE);

  fprintf (stderr, "ERROR:  Hash table full\n");
  assert (FALSE);
}




//  Insert string subscript  i  into the global hash table.
//  Sequence and information about the string are in
//  global variables  basesData, String_Start, String_Info, ....
static
void
Put_String_In_Hash(uint32 curID, uint32 i) {
  String_Ref_t  ref = 0;
  int           skip_ct;
  uint64        key;
  uint64        key_is_bad;
  int           j;

  uint32        kmers_skipped  = 0;
  uint32        kmers_bad      = 0;
  uint32        kmers_inserted = 0;

  char *p      = basesData + String_Start[i];
  char *window = basesData + String_Start[i];

  key = key_is_bad = 0;

  for (uint32 j=0;  j<G.Kmer_Len; j ++) {
    key_is_bad |= (uint64) (Char_Is_Bad[(int) * p]) << j;
    key        |= (uint64) (Bit_Equivalent[(int) * (p ++)]) << (2 * j);
  }

  setStringRefStringNum(ref, i);

  if (i > MAX_STRING_NUM)
    fprintf (stderr, "Too many strings for hash table--exiting\n"), exit(1);

  setStringRefOffset(ref, TRUELY_ZERO);

  skip_ct = 0;

  setStringRefEmpty(ref, TRUELY_ZERO);

  if (key_is_bad == false) {
    Hash_Insert(ref, key, window);
    kmers_inserted++;

  } else {
    kmers_bad++;
  }

  while (*p != 0) {
    window++;

    String_Ref_t newoff = getStringRefOffset(ref) + 1;
    assert(newoff < OFFSET_MASK);

    setStringRefOffset(ref, newoff);

    if (++skip_ct > HASH_KMER_SKIP)
      skip_ct = 0;

    key_is_bad >>= 1;
    key_is_bad |= (uint64) (Char_Is_Bad[(int) * p]) << (G.Kmer_Len - 1);

    key >>= 2;
    key  |= (uint64) (Bit_Equivalent[(int) * (p ++)]) << (2 * (G.Kmer_Len - 1));

    if (skip_ct > 0) {
      kmers_skipped++;
      continue;
    }

    if (key_is_bad) {
      kmers_bad++;
      continue;
    }

    Hash_Insert(ref, key, window);
    kmers_inserted++;
  }

  //fprintf(stderr, "STRING %u skipped %u bad %u inserted %u\n",
  //        curID, kmers_skipped, kmers_bad, kmers_inserted);
}



// Read the next batch of strings from  stream  and create a hash
//  table index of their  G.Kmer_Len -mers.  Return  1  if successful;
//  0 otherwise.  The batch ends when either end-of-file is encountered
//  or  Max_Hash_Strings  have been read in.   first_frag_id  is the
//  internal ID of the first fragment in the hash table.
int
Build_Hash_Index(gkStore *gkpStore, uint32 bgnID, uint32 endID) {
  String_Ref_t  ref;
  uint64  total_len;
  uint64   hash_entry_limit;

  fprintf(stderr, "Build_Hash_Index from "F_U32" to "F_U32"\n", bgnID, endID);

  Hash_String_Num_Offset = bgnID;
  String_Ct              = 0;
  Extra_String_Ct        = 0;
  Extra_String_Subcount  = MAX_EXTRA_SUBCOUNT;
  total_len              = 0;

  //if (Data == NULL) {
  //  Extra_Data_Len    = Max_Hash_Data_Len + AS_MAX_READLEN;
  //  Data_Len          = Max_Hash_Data_Len + AS_MAX_READLEN;
  //
  //  basesData         = new char [Data_Len];
  //  qualsData         = new char [Data_Len];
  //
  //  old_ref_len       = Data_Len / (HASH_KMER_SKIP + 1);
  //  nextRef           = new String_Ref_t [old_ref_len];
  //}

  //memset(nextRef,         0xff, old_ref_len     * sizeof(String_Ref_t));

  memset(Hash_Table,       0x00, HASH_TABLE_SIZE * sizeof(Hash_Bucket_t));
  memset(Hash_Check_Array, 0x00, HASH_TABLE_SIZE * sizeof(Check_Vector_t));

  Extra_Ref_Ct     = 0;
  Hash_Entries     = 0;
  hash_entry_limit = G.Max_Hash_Load * HASH_TABLE_SIZE * ENTRIES_PER_BUCKET;

#if 0
  fprintf(stderr, "HASH LOADING STARTED: fragID   %12"F_U64P"\n", first_frag_id);
  fprintf(stderr, "HASH LOADING STARTED: strings  %12"F_U64P" out of %12"F_U64P" max.\n", String_Ct, G.Max_Hash_Strings);
  fprintf(stderr, "HASH LOADING STARTED: length   %12"F_U64P" out of %12"F_U64P" max.\n", total_len, G.Max_Hash_Data_Len);
  fprintf(stderr, "HASH LOADING STARTED: entries  %12"F_U64P" out of %12"F_U64P" max (load %.2f).\n", Hash_Entries, hash_entry_limit,
         (100.0 * Hash_Entries) / (HASH_TABLE_SIZE * ENTRIES_PER_BUCKET));
#endif

  //  Compute an upper limit on the number of bases we will load.  The number of Hash_Entries
  //  can't be computed here, so the real loop below could end earlier than expected - and we
  //  don't use a little bit of memory.

  uint32  nSkipped  = 0;
  uint32  nShort    = 0;
  uint32  nLoadable = 0;

  uint64  maxAlloc = 0;
  uint32  curID    = 0;  //  The last ID loaded into the hash

  for (curID=bgnID; ((String_Ct <  G.Max_Hash_Strings) &&
                     (total_len <  G.Max_Hash_Data_Len) &&
                     (curID     <= endID)); curID++) {
    gkRead *read = gkpStore->gkStore_getRead(curID);

    if ((read->gkRead_libraryID() < G.minLibToHash) ||
        (read->gkRead_libraryID() > G.maxLibToHash)) {
      nSkipped++;
      continue;
    }

    if (read->gkRead_sequenceLength() < G.Min_Olap_Len) {
      nShort++;
      continue;
    }

    nLoadable++;

    maxAlloc += read->gkRead_sequenceLength() + 1;
  }

  fprintf(stderr, "Found "F_U32" reads with length "F_U64" to load; "F_U32" skipped by being too short; "F_U32" skipped per library restriction\n",
          nLoadable, maxAlloc, nShort, nSkipped);

  //  This should be less than what the user requested on the command line

  if (maxAlloc >= G.Max_Hash_Data_Len + AS_MAX_READLEN)
    fprintf(stderr, "maxAlloc = "F_U64" G.Max_Hash_Data_Len = "F_U64"  AS_MAX_READLEN = %u\n", maxAlloc, G.Max_Hash_Data_Len, AS_MAX_READLEN);
  assert(maxAlloc < G.Max_Hash_Data_Len + AS_MAX_READLEN);

  //  Allocate space, then fill it.

  uint64 nextRef_Len = maxAlloc / (HASH_KMER_SKIP + 1);
  Extra_Data_Len = Data_Len  = maxAlloc;

  basesData = new char         [Data_Len];
  qualsData = new char         [Data_Len];
  nextRef   = new String_Ref_t [nextRef_Len];

  memset(nextRef, 0xff, sizeof(String_Ref_t) * nextRef_Len);

  gkReadData   *readData = new gkReadData;

  for (curID=bgnID; ((String_Ct    <  G.Max_Hash_Strings) &&
                     (total_len    <  G.Max_Hash_Data_Len) &&
                     (Hash_Entries <  hash_entry_limit) &&
                     (curID        <= endID)); curID++, String_Ct++) {

    //  Load sequence if it exists, otherwise, add an empty read.
    //  Duplicated in Process_Overlaps().

    String_Start[String_Ct]                    = UINT64_MAX;

    String_Info[String_Ct].length              = 0;
    String_Info[String_Ct].lfrag_end_screened  = TRUE;
    String_Info[String_Ct].rfrag_end_screened  = TRUE;

    gkRead  *read = gkpStore->gkStore_getRead(curID);

    if ((read->gkRead_libraryID() < G.minLibToHash) ||
        (read->gkRead_libraryID() > G.maxLibToHash))
      continue;

    uint32 len = read->gkRead_sequenceLength();

    if (len < G.Min_Olap_Len)
      continue;

    gkpStore->gkStore_loadReadData(read, readData);

    char   *seqptr   = readData->gkReadData_getSequence();
    char   *qltptr   = readData->gkReadData_getQualities();

    //  Note where we are going to store the string, and how long it is

    String_Start[String_Ct]                    = total_len;

    String_Info[String_Ct].length              = len;
    String_Info[String_Ct].lfrag_end_screened  = FALSE;
    String_Info[String_Ct].rfrag_end_screened  = FALSE;

    //  Store it.

    for (uint32 i=0; i<len; i++, total_len++) {
      basesData[total_len] = tolower(seqptr[i]);
      qualsData[total_len] = qltptr[i];
    }

    basesData[total_len] = 0;
    qualsData[total_len] = 0;

    total_len++;

    //  Skipping kners is totally untested.
#if 0
    if (HASH_KMER_SKIP > 0) {
      uint32 extra   = new_len % (HASH_KMER_SKIP + 1);

      if (extra > 0)
        new_len += 1 + HASH_KMER_SKIP - extra;
    }
#endif

    //  Trouble - allocate more space for sequence and quality data.
    //  This was computed ahead of time!

    if (total_len > maxAlloc)
      fprintf(stderr, "total_len="F_U64"  len="F_U32"  maxAlloc="F_U64"\n", total_len, len, maxAlloc);
    assert(total_len <= maxAlloc);

    //  What is Extra_Data_Len?  It's set to Data_Len if we would have reallocated here.

    Put_String_In_Hash(curID, String_Ct);

    if ((String_Ct % 100000) == 0)
      fprintf (stderr, "String_Ct:%12"F_U64P"/%12"F_U32P"  totalLen:%12"F_U64P"/%12"F_U64P"  Hash_Entries:%12"F_U64P"/%12"F_U64P"  Load: %.2f%%\n",
               String_Ct,    G.Max_Hash_Strings,
               total_len,    G.Max_Hash_Data_Len,
               Hash_Entries,
               hash_entry_limit,
               100.0 * Hash_Entries / (HASH_TABLE_SIZE * ENTRIES_PER_BUCKET));
  }

  curID--;  //  We always stop on the read after we loaded.

  delete readData;

  fprintf(stderr, "HASH LOADING STOPPED: strings  %12"F_U64P" out of %12"F_U32P" max.\n", String_Ct, G.Max_Hash_Strings);
  fprintf(stderr, "HASH LOADING STOPPED: length   %12"F_U64P" out of %12"F_U64P" max.\n", total_len, G.Max_Hash_Data_Len);
  fprintf(stderr, "HASH LOADING STOPPED: entries  %12"F_U64P" out of %12"F_U64P" max (load %.2f).\n", Hash_Entries, hash_entry_limit,
          100.0 * Hash_Entries / (HASH_TABLE_SIZE * ENTRIES_PER_BUCKET));

  if (String_Ct == 0) {
    fprintf(stderr, "HASH LOADING STOPPED: no strings added?\n");
    return(endID);
  }

  Used_Data_Len = total_len;

  //fprintf(stderr, "Extra_Ref_Ct = "F_U64"  Max_Extra_Ref_Space = "F_U64"\n", Extra_Ref_Ct, Max_Extra_Ref_Space);

  if (Extra_Ref_Ct > Max_Extra_Ref_Space) {
    int32          newSize  = (Max_Extra_Ref_Space == 0) ? 16 * 1024 : Max_Extra_Ref_Space * 2;

    while (newSize < Extra_Ref_Ct)
      newSize *= 2;

    String_Ref_t  *newSpace = new String_Ref_t [newSize];

    memcpy(newSpace, Extra_Ref_Space, sizeof(String_Ref_t) * Max_Extra_Ref_Space);

    delete [] Extra_Ref_Space;

    Max_Extra_Ref_Space = newSize;    //  Former max_extra_ref_ct
    Extra_Ref_Space     = newSpace;
  }


  if (G.Kmer_Skip_File != NULL)
    Mark_Skip_Kmers();


  // Coalesce reference chain into adjacent entries in  Extra_Ref_Space
  Extra_Ref_Ct = 0;
  for (int32 i = 0;  i < HASH_TABLE_SIZE;  i ++)
    for (int32 j = 0;  j < Hash_Table[i].Entry_Ct;  j ++) {
      ref = Hash_Table[i].Entry[j];
      if (! getStringRefLast(ref) && ! getStringRefEmpty(ref)) {
        Extra_Ref_Space[Extra_Ref_Ct] = ref;
        setStringRefStringNum(Hash_Table[i].Entry[j], (String_Ref_t)(Extra_Ref_Ct >> OFFSET_BITS));
        setStringRefOffset  (Hash_Table[i].Entry[j], (String_Ref_t)(Extra_Ref_Ct & OFFSET_MASK));
        Extra_Ref_Ct ++;
        do {
          ref = nextRef[(String_Start[getStringRefStringNum(ref)] + getStringRefOffset(ref)) / (HASH_KMER_SKIP + 1)];
          Extra_Ref_Space[Extra_Ref_Ct ++] = ref;
        }  while (! getStringRefLast(ref));
      }
    }

  return(curID);
}
