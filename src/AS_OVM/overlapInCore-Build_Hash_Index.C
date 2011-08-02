
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

static const char *rcsid = "$Id: overlapInCore-Build_Hash_Index.C,v 1.2 2011-08-02 02:21:03 brianwalenz Exp $";

#include "overlapInCore.H"




//  Add string  s  as an extra hash table string and return
//  a single reference to the beginning of it.
static
String_Ref_t
Add_Extra_Hash_String(const char * s) {
  String_Ref_t  ref = 0;
  String_Ref_t  sub = 0;
  size_t  new_len;
  int  len;
   
  new_len = Used_Data_Len + Kmer_Len;
  if  (Extra_String_Subcount < MAX_EXTRA_SUBCOUNT) {
    sub = String_Ct + Extra_String_Ct - 1;
  } else {
    sub = String_Ct + Extra_String_Ct;
    if  (sub >= String_Start_Size) {
      String_Start_Size *= MEMORY_EXPANSION_FACTOR;
      if  (sub >= String_Start_Size)
        String_Start_Size = sub;
      String_Start = (int64 *) safe_realloc (String_Start,
                                             String_Start_Size * sizeof (int64));
    }
    String_Start [sub] = Used_Data_Len;
    Extra_String_Ct ++;
    Extra_String_Subcount = 0;
    new_len ++;
  }

  if  (new_len >= Extra_Data_Len) {
    Extra_Data_Len = (size_t) (Extra_Data_Len * MEMORY_EXPANSION_FACTOR);
    if  (new_len > Extra_Data_Len)
      Extra_Data_Len = new_len;
    Data = (char *) safe_realloc (Data, Extra_Data_Len);
  }
  strncpy (Data + String_Start [sub] + Kmer_Len * Extra_String_Subcount,
           s, Kmer_Len + 1);
  Used_Data_Len = new_len;

  setStringRefStringNum(ref, sub);
  if  (sub > MAX_STRING_NUM) {
    fprintf (stderr, "Too many skip kmer strings for hash table.\n"
             "Try skipping hopeless check (-z option)\n"
             "Exiting\n");
    exit (1);
  }
  setStringRefOffset(ref, (String_Ref_t)Extra_String_Subcount * (String_Ref_t)Kmer_Len);

  assert(Extra_String_Subcount * Kmer_Len < OFFSET_MASK);

  setStringRefLast(ref, TRUELY_ONE);
  setStringRefEmpty(ref, TRUELY_ONE);

  Extra_String_Subcount ++;

  return  ref;
}





//  Mark  left/right_end_screened in global  String_Info  for
//   ref  and everything in its list, if they occur near
//  enough to the end of the string.

static
void
Mark_Screened_Ends_Single(String_Ref_t ref) {

  int32 s_num = getStringRefStringNum(ref);
  int32 len = String_Info [s_num] . length;

  if  (getStringRefOffset(ref) < HOPELESS_MATCH) {
#ifdef COMPARE
    fprintf(stderr, "Mark_Screened_Ends_Single()-- %d mark left\n", s_num);
#endif
    String_Info [s_num] . lfrag_end_screened = TRUE;
  }

  if  (len - getStringRefOffset(ref) - Kmer_Len + 1 < HOPELESS_MATCH) {
#ifdef COMPARE
    fprintf(stderr, "Mark_Screened_Ends_Single()-- %d mark left\n", s_num);
#endif
    String_Info [s_num] . rfrag_end_screened = TRUE;
  }
}



static
void
Mark_Screened_Ends_Chain(String_Ref_t ref) {

  Mark_Screened_Ends_Single (ref);

  while  (! getStringRefLast(ref)) {
    ref = Next_Ref [(String_Start [getStringRefStringNum(ref)] + getStringRefOffset(ref)) / (HASH_KMER_SKIP + 1)];
    Mark_Screened_Ends_Single (ref);
  }
}


//  Set the  empty  bit to true for the hash table entry
//  corresponding to string  s  whose hash key is  key .
//  Also set global  String_Info . left/right_end_screened
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
    for  (i = 0;  i < Hash_Table [sub] . Entry_Ct;  i ++)
      if  (Hash_Table [sub] . Check [i] == key_check) {
        h_ref = Hash_Table [sub] . Entry [i];
        t = Data + String_Start [getStringRefStringNum(h_ref)] + getStringRefOffset(h_ref);
        if  (strncmp (s, t, Kmer_Len) == 0) {
          if  (! getStringRefEmpty(Hash_Table [sub] . Entry [i]))
            Mark_Screened_Ends_Chain (Hash_Table [sub] . Entry [i]);
          setStringRefEmpty(Hash_Table [sub] . Entry [i], TRUELY_ONE);
          return;
        }
      }
    assert (i == Hash_Table [sub] . Entry_Ct);
    if  (Hash_Table [sub] . Entry_Ct < ENTRIES_PER_BUCKET) {
      // Not found
      if  (Use_Hopeless_Check) {
        Hash_Table [sub] . Entry [i] = Add_Extra_Hash_String (s);
        setStringRefEmpty(Hash_Table [sub] . Entry [i], TRUELY_ONE);
        Hash_Table [sub] . Check [i] = key_check;
        Hash_Table [sub] . Entry_Ct ++;
        Hash_Table [sub] . Hits [i] = 0;
        Hash_Entries ++;
        shift = HASH_CHECK_FUNCTION (key);
        Hash_Check_Array [sub] |= (((Check_Vector_t) 1) << shift);
      }
      return;
    }
    sub = (sub + probe) % HASH_TABLE_SIZE;
  }  while  (++ ct < HASH_TABLE_SIZE);

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
  char  line [MAX_LINE_LEN];
  int  ct = 0;

  rewind (Kmer_Skip_File);

  while  (fgets (line, MAX_LINE_LEN, Kmer_Skip_File) != NULL) {
    int  i, len;

    ct ++;
    len = strlen (line) - 1;
    if  (line [0] != '>' || line [len] != '\n') {
      fprintf (stderr, "ERROR:  Bad line %d in kmer skip file\n", ct);
      fputs (line, stderr);
      exit (1);
    }

    if  (fgets (line, MAX_LINE_LEN, Kmer_Skip_File) == NULL) {
      fprintf (stderr, "ERROR:  Bad line after %d in kmer skip file\n", ct);
      exit (1);
    }
    ct ++;
    len = strlen (line) - 1;
    if  (len != Kmer_Len || line [len] != '\n') {
      fprintf (stderr, "ERROR:  Bad line %d in kmer skip file\n", ct);
      fputs (line, stderr);
      exit (1);
    }
    line [len] = '\0';

    key = 0;
    for  (i = 0;  i < len;  i ++) {
      line [i] = tolower (line [i]);
      key |= (uint64) (Bit_Equivalent [(int) line [i]]) << (2 * i);
    }
    Hash_Mark_Empty (key, line);

    reverseComplementSequence (line, len);
    key = 0;
    for  (i = 0;  i < len;  i ++)
      key |= (uint64) (Bit_Equivalent [(int) line [i]]) << (2 * i);
    Hash_Mark_Empty (key, line);
  }

  fprintf (stderr, "String_Ct = %d  Extra_String_Ct = %d  Extra_String_Subcount = %d\n",
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
  Hash_Check_Array [Sub] |= (((Check_Vector_t) 1) << Shift);
  Key_Check = KEY_CHECK_FUNCTION (Key);
  Probe = PROBE_FUNCTION (Key);

  Ct = 0;
  do {
    for  (i = 0;  i < Hash_Table [Sub] . Entry_Ct;  i ++)
      if  (Hash_Table [Sub] . Check [i] == Key_Check) {
        H_Ref = Hash_Table [Sub] . Entry [i];
        T = Data + String_Start [getStringRefStringNum(H_Ref)] + getStringRefOffset(H_Ref);
        if  (strncmp (S, T, Kmer_Len) == 0) {
          if  (getStringRefLast(H_Ref)) {
            Extra_Ref_Ct ++;
          }
          Next_Ref [(String_Start [getStringRefStringNum(Ref)] + getStringRefOffset(Ref)) / (HASH_KMER_SKIP + 1)] = H_Ref;
          Extra_Ref_Ct ++;
          setStringRefLast(Ref, TRUELY_ZERO);
          Hash_Table [Sub] . Entry [i] = Ref;

          if  (Hash_Table [Sub] . Hits [i] < HIGHEST_KMER_LIMIT)
            Hash_Table [Sub] . Hits [i] ++;

          return;
        }
      }
    if  (i != Hash_Table [Sub] . Entry_Ct) {
      fprintf (stderr, "i = %d  Sub = " F_S64 "  Entry_Ct = %d\n",
               i, Sub, Hash_Table [Sub] . Entry_Ct);
    }
    assert (i == Hash_Table [Sub] . Entry_Ct);
    if  (Hash_Table [Sub] . Entry_Ct < ENTRIES_PER_BUCKET) {
      setStringRefLast(Ref, TRUELY_ONE);
      Hash_Table [Sub] . Entry [i] = Ref;
      Hash_Table [Sub] . Check [i] = Key_Check;
      Hash_Table [Sub] . Entry_Ct ++;
      Hash_Entries ++;
      Hash_Table [Sub] . Hits [i] = 1;
      return;
    }
    Sub = (Sub + Probe) % HASH_TABLE_SIZE;
  }  while  (++ Ct < HASH_TABLE_SIZE);

  fprintf (stderr, "ERROR:  Hash table full\n");
  assert (FALSE);
}




//  Insert string subscript  i  into the global hash table.
//  Sequence and information about the string are in
//  global variables  Data, String_Start, String_Info, ....
static
void
Put_String_In_Hash(int i) {
  String_Ref_t  ref = 0;
  char  * p, * window;
  int  kmers_inserted = 0;
  int  skip_ct;
  uint64  key, key_is_bad;
  int  j;

  if  (String_Info [i] . length < Kmer_Len)
    return;

  p = window = Data + String_Start [i];
  key = key_is_bad = 0;
  for  (j = 0;  j < Kmer_Len;  j ++) {
    key_is_bad |= (uint64) (Char_Is_Bad [(int) * p]) << j;
    key |= (uint64) (Bit_Equivalent [(int) * (p ++)]) << (2 * j);
  }

  setStringRefStringNum(ref, i);
  if  (i > MAX_STRING_NUM) {
    fprintf (stderr, "Too many strings for hash table--exiting\n");
    exit (1);
  }
  setStringRefOffset(ref, TRUELY_ZERO);
  skip_ct = 0;
  setStringRefEmpty(ref, TRUELY_ZERO);

  if  (key_is_bad == false) {
    Hash_Insert (ref, key, window);
    kmers_inserted ++;
  }

  while  ((* p) != '\0') {
    window ++;

    {
      String_Ref_t newoff = getStringRefOffset(ref) + 1;
      assert(newoff < OFFSET_MASK);
      setStringRefOffset(ref, newoff);
    }

    if  (++ skip_ct > HASH_KMER_SKIP)
      skip_ct = 0;

    key_is_bad >>= 1;
    key_is_bad |= (uint64) (Char_Is_Bad [(int) * p]) << (Kmer_Len - 1);
    key >>= 2;
    key |= (uint64) (Bit_Equivalent [(int) * (p ++)]) << (2 * (Kmer_Len - 1));

    if  (skip_ct == 0 && ! key_is_bad) {
      Hash_Insert (ref, key, window);
      kmers_inserted ++;
    }
  }
}



// Read the next batch of strings from  stream  and create a hash
//  table index of their  Kmer_Len -mers.  Return  1  if successful;
//  0 otherwise.  The batch ends when either end-of-file is encountered
//  or  Max_Hash_Strings  have been read in.   first_frag_id  is the
//  internal ID of the first fragment in the hash table.
int
Build_Hash_Index(gkStream *stream, int32 first_frag_id, gkFragment *myRead) {
  String_Ref_t  ref;
  int64  total_len;
  static int64  max_extra_ref_ct = 0;
  static int64  old_ref_len, new_ref_len;
  int  frag_status;
  int64  i;
  int  hash_entry_limit;
  int  j;

  Hash_String_Num_Offset = first_frag_id;
  String_Ct = Extra_String_Ct = 0;
  Extra_String_Subcount = MAX_EXTRA_SUBCOUNT;
  total_len = 0;
  if  (Data == NULL) {
    Extra_Data_Len = Data_Len = Max_Hash_Data_Len + AS_READ_MAX_NORMAL_LEN;
    Data = (char *) safe_realloc (Data, Data_Len);
    Quality_Data = (char *) safe_realloc (Quality_Data, Data_Len);
    old_ref_len = Data_Len / (HASH_KMER_SKIP + 1);
    Next_Ref = (String_Ref_t *) safe_realloc
      (Next_Ref, old_ref_len * sizeof (String_Ref_t));
  }

  memset (Next_Ref, '\377', old_ref_len * sizeof (String_Ref_t));
  memset (Hash_Table, 0, HASH_TABLE_SIZE * sizeof (Hash_Bucket_t));
  memset (Hash_Check_Array, 0, HASH_TABLE_SIZE * sizeof (Check_Vector_t));

  fprintf (stderr, "### Build_Hash:  first_frag_id = %d  Max_Hash_Strings = %d\n",
           first_frag_id, Max_Hash_Strings);

  Extra_Ref_Ct = 0;
  Hash_Entries = 0;
  hash_entry_limit = Max_Hash_Load * HASH_TABLE_SIZE * ENTRIES_PER_BUCKET;

  while  (String_Ct < Max_Hash_Strings
          && total_len < Max_Hash_Data_Len
          && Hash_Entries < hash_entry_limit
          && (frag_status
              = Read_Next_Frag (Sequence_Buffer, Quality_Buffer, stream,
                                myRead, & Last_Hash_Frag_Read, minLibToHash, maxLibToHash))) {
    int  extra, len;
    size_t  new_len;

    if  (frag_status == DELETED_FRAG) {
      Sequence_Buffer [0] = '\0';
      Quality_Buffer [0] = '\0';
    }

    String_Start [String_Ct] = total_len;
    len = strlen (Sequence_Buffer);
    String_Info [String_Ct] . length = len;
    String_Info [String_Ct] . lfrag_end_screened = FALSE;
    String_Info [String_Ct] . rfrag_end_screened = FALSE;
    new_len = total_len + len + 1;
    extra = new_len % (HASH_KMER_SKIP + 1);
    if  (extra > 0)
      new_len += 1 + HASH_KMER_SKIP - extra;

    if  (new_len > Data_Len) {
      Data_Len = (size_t) (Data_Len * MEMORY_EXPANSION_FACTOR);
      if  (new_len > Data_Len)
        Data_Len = new_len;
      if  (Data_Len > Extra_Data_Len) {
        Data = (char *) safe_realloc (Data, Data_Len);
        Extra_Data_Len = Data_Len;
      }
      Quality_Data = (char *) safe_realloc (Quality_Data, Data_Len);
      new_ref_len = Data_Len / (HASH_KMER_SKIP + 1);
      Next_Ref = (String_Ref_t *) safe_realloc
        (Next_Ref, new_ref_len * sizeof (String_Ref_t));
      memset (Next_Ref + old_ref_len, '\377',
              (new_ref_len - old_ref_len) * sizeof (String_Ref_t));
      old_ref_len = new_ref_len;
    }

    strcpy (Data + total_len, Sequence_Buffer);
    memcpy (Quality_Data + total_len, Quality_Buffer, len + 1);
    total_len = new_len;

    Put_String_In_Hash (String_Ct);

    if ((String_Ct % 100000) == 0)
      fprintf (stderr, "String_Ct:%d  totalLen:"F_S64"  Hash_Entries:"F_S64"  Load:%.1f%%\n",
               String_Ct,
               total_len,
               Hash_Entries,
               (100.0 * Hash_Entries) / (HASH_TABLE_SIZE * ENTRIES_PER_BUCKET));

    String_Ct ++;
  }

  if  (String_Ct == 0)
    return  0;

  fprintf (stderr, "strings read = %d  total_len = " F_S64 "\n",
           String_Ct, total_len);
  Used_Data_Len = total_len;

  fprintf (stderr, "Hash_Entries = " F_S64 "  Load = %.1f%%\n",
           Hash_Entries, (100.0 * Hash_Entries) / (HASH_TABLE_SIZE * ENTRIES_PER_BUCKET));

  if  (Extra_Ref_Ct > max_extra_ref_ct) {
    max_extra_ref_ct *= MEMORY_EXPANSION_FACTOR;
    if  (Extra_Ref_Ct > max_extra_ref_ct)
      max_extra_ref_ct = Extra_Ref_Ct;
    fprintf (stderr,
             "### realloc  Extra_Ref_Space  max_extra_ref_ct = " F_S64 "\n",
             max_extra_ref_ct);
    Extra_Ref_Space = (String_Ref_t *) safe_realloc (Extra_Ref_Space,
                                                     max_extra_ref_ct * sizeof (String_Ref_t));
  }

  if  (Kmer_Skip_File != NULL)
    Mark_Skip_Kmers ();

  // Coalesce reference chain into adjacent entries in  Extra_Ref_Space
  Extra_Ref_Ct = 0;
  for  (i = 0;  i < HASH_TABLE_SIZE;  i ++)
    for  (j = 0;  j < Hash_Table [i] . Entry_Ct;  j ++) {
      ref = Hash_Table [i] . Entry [j];
      if  (! getStringRefLast(ref) && ! getStringRefEmpty(ref)) {
        Extra_Ref_Space [Extra_Ref_Ct] = ref;
        setStringRefStringNum(Hash_Table [i] . Entry [j], (String_Ref_t)(Extra_Ref_Ct >> OFFSET_BITS));
        setStringRefOffset   (Hash_Table [i] . Entry [j], (String_Ref_t)(Extra_Ref_Ct & OFFSET_MASK));
        Extra_Ref_Ct ++;
        do {
          ref = Next_Ref [(String_Start [getStringRefStringNum(ref)] + getStringRefOffset(ref)) / (HASH_KMER_SKIP + 1)];
          Extra_Ref_Space [Extra_Ref_Ct ++] = ref;
        }  while  (! getStringRefLast(ref));
      }
    }

  return  1;
}
