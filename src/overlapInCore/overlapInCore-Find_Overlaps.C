
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
Add_Match(String_Ref_t ref,
          int * start,
          int offset,
          int * consistent,
          Work_Area_t * WA) {
  int  * p, save;
  int  diag = 0, new_diag, expected_start = 0, num_checked = 0;
  int  move_to_front = false;

  new_diag = getStringRefOffset(ref) - offset;

  for (p = start;  (* p) != 0;  p = & (WA->Match_Node_Space [(* p)].Next)) {
    expected_start = WA->Match_Node_Space [(* p)].Start + WA->Match_Node_Space [(* p)].Len - G.Kmer_Len + 1 + HASH_KMER_SKIP;

    diag = WA->Match_Node_Space [(* p)].Offset - WA->Match_Node_Space [(* p)].Start;

    if (expected_start < offset)
      break;

    if (expected_start == offset) {
      if (new_diag == diag) {
        WA->Match_Node_Space [(* p)].Len += 1 + HASH_KMER_SKIP;
        if (move_to_front) {
          save = (* p);
          (* p) = WA->Match_Node_Space [(* p)].Next;
          WA->Match_Node_Space [save].Next = (* start);
          (* start) = save;
        }
        return;
      } else
        move_to_front = true;
    }
    num_checked ++;
  }

  if (WA->Next_Avail_Match_Node == WA->Match_Node_Size) {
    int32          newSize  = WA->Match_Node_Size * 2;
    Match_Node_t  *newSpace = new Match_Node_t [newSize];

    memcpy(newSpace, WA->Match_Node_Space, sizeof(Match_Node_t) * WA->Match_Node_Size);

    delete [] WA->Match_Node_Space;

    WA->Match_Node_Size  = newSize;
    WA->Match_Node_Space = newSpace;
  }

  if ((* start) != 0
      && (num_checked > 0
          || abs (diag - new_diag) > 3
          || offset < expected_start + G.Kmer_Len - 2))
    (* consistent) = false;

  save = (* start);
  (* start) = WA->Next_Avail_Match_Node;
  WA->Next_Avail_Match_Node ++;

  WA->Match_Node_Space [(* start)].Offset = getStringRefOffset(ref);
  WA->Match_Node_Space [(* start)].Len = G.Kmer_Len;
  WA->Match_Node_Space [(* start)].Start = offset;
  WA->Match_Node_Space [(* start)].Next = save;

#if 0
  fprintf(stderr, "Add_Match()-- %3d offset %d len %d start %d next %d\n",
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
static
void
Add_Ref(String_Ref_t Ref, int Offset, Work_Area_t * WA) {
  uint32  Prev, StrNum, Sub;
  int  consistent;

  StrNum = getStringRefStringNum(Ref);
  Sub = (StrNum ^ (StrNum >> STRING_OLAP_SHIFT)) & STRING_OLAP_MASK;

  while (WA->String_Olap_Space [Sub].Full
          && WA->String_Olap_Space [Sub].String_Num != StrNum) {
    Prev = Sub;
    Sub = WA->String_Olap_Space [Sub].Next;
    if (Sub == 0) {
      if (WA->Next_Avail_String_Olap == WA->String_Olap_Size) {
        int32          newSize  = WA->String_Olap_Size * 2;
        String_Olap_t *newSpace = new String_Olap_t [newSize];

        memcpy(newSpace, WA->String_Olap_Space, sizeof(String_Olap_t) * WA->String_Olap_Size);

        delete [] WA->String_Olap_Space;

        WA->String_Olap_Size  = newSize;
        WA->String_Olap_Space = newSpace;
      }

      Sub = WA->Next_Avail_String_Olap ++;
      WA->String_Olap_Space [Prev].Next = Sub;
      WA->String_Olap_Space [Sub].Full = false;
      break;
    }
  }

  if (! WA->String_Olap_Space [Sub].Full) {
    WA->String_Olap_Space [Sub].String_Num = StrNum;
    WA->String_Olap_Space [Sub].Match_List = 0;
    WA->String_Olap_Space [Sub].diag_sum = 0.0;
    WA->String_Olap_Space [Sub].diag_ct = 0;
    WA->String_Olap_Space [Sub].diag_bgn = AS_MAX_READLEN;
    WA->String_Olap_Space [Sub].diag_end = 0;
    WA->String_Olap_Space [Sub].Next = 0;
    WA->String_Olap_Space [Sub].Full = true;
    WA->String_Olap_Space [Sub].consistent = true;
  }

  consistent = WA->String_Olap_Space [Sub].consistent;

  WA->String_Olap_Space [Sub].diag_sum += (double)getStringRefOffset(Ref) - Offset;
  WA->String_Olap_Space [Sub].diag_ct ++;
  if (WA->String_Olap_Space [Sub].diag_bgn > Offset) WA->String_Olap_Space [Sub].diag_bgn = Offset;
  if (WA->String_Olap_Space [Sub].diag_end < Offset) WA->String_Olap_Space [Sub].diag_end = Offset;
  Add_Match (Ref, & (WA->String_Olap_Space [Sub].Match_List), Offset, & consistent, WA);

  WA->String_Olap_Space [Sub].consistent = consistent;

  return;
}




//  Search for string  S  with hash key  Key  in the global
//  Hash_Table  starting at subscript  Sub. Return the matching
//  reference in the hash table if there is one, or else a reference
//  with the  Empty bit set true.  Set  (* Where)  to the subscript in
//  Extra_Ref_Space  where the reference was found if it was found there.
//  Set  (* hi_hits)  to  true  if hash table entry is found but is empty
//  because it was screened out, otherwise set to false.
static
String_Ref_t
Hash_Find(uint64 Key, int64 Sub, char * S, int64 * Where, int * hi_hits) {
  String_Ref_t  H_Ref = 0;
  char  * T;
  unsigned char  Key_Check;
  int64  Ct, Probe;
  int  i;

  Key_Check = KEY_CHECK_FUNCTION (Key);
  Probe = PROBE_FUNCTION (Key);

  (* hi_hits) = false;
  Ct = 0;
  do {
    for (i = 0;  i < Hash_Table [Sub].Entry_Ct;  i ++)
      if (Hash_Table [Sub].Check [i] == Key_Check) {
        int  is_empty;

        H_Ref = Hash_Table [Sub].Entry [i];
        //fprintf(stderr, "Href = Hash_Table %u Entry %u = " F_U64 "\n", Sub, i, H_Ref);

        is_empty = getStringRefEmpty(H_Ref);
        if (! getStringRefLast(H_Ref) && ! is_empty) {
          (* Where) = ((uint64)getStringRefStringNum(H_Ref) << OFFSET_BITS) + getStringRefOffset(H_Ref);
          H_Ref = Extra_Ref_Space [(* Where)];
          //fprintf(stderr, "Href = Extra_Ref_Space " F_U64 " = " F_U64 "\n", *Where, H_Ref);
        }
        //fprintf(stderr, "Href = " F_U64 "  Get String_Start[ " F_U64 " ] + " F_U64 "\n", getStringRefStringNum(H_Ref), getStringRefOffset(H_Ref));
        T = basesData + String_Start [getStringRefStringNum(H_Ref)] + getStringRefOffset(H_Ref);
        if (strncmp (S, T, G.Kmer_Len) == 0) {
          if (is_empty) {
            setStringRefEmpty(H_Ref, TRUELY_ONE);
            (* hi_hits) = true;
          }
          return  H_Ref;
        }
      }
    if (Hash_Table [Sub].Entry_Ct < ENTRIES_PER_BUCKET) {
      setStringRefEmpty(H_Ref, TRUELY_ONE);
      return  H_Ref;
    }
    Sub = (Sub + Probe) % HASH_TABLE_SIZE;
  }  while (++ Ct < HASH_TABLE_SIZE);

  setStringRefEmpty(H_Ref, TRUELY_ONE);
  return  H_Ref;
}






//  Find and output all overlaps and branch points between string
//   Frag  and any fragment currently in the global hash table.
//   Frag_Len  is the length of  Frag  and  Frag_Num  is its ID number.
//   Dir  is the orientation of  Frag .

void
Find_Overlaps(char Frag [], int Frag_Len, uint32 Frag_Num, Direction_t Dir, Work_Area_t * WA) {
  String_Ref_t  Ref;
  char  * P, * Window;
  uint64  Key, Next_Key;
  int64  Sub, Next_Sub, Where;
  Check_Vector_t  This_Check, Next_Check;
  int  Offset, Shift, Next_Shift;
  int  hi_hits;
  int  j;

  memset (WA->String_Olap_Space, 0, STRING_OLAP_MODULUS * sizeof (String_Olap_t));
  WA->Next_Avail_String_Olap = STRING_OLAP_MODULUS;
  WA->Next_Avail_Match_Node = 1;

  assert (Frag_Len >= G.Kmer_Len);

  Offset = 0;
  P = Window = Frag;

  WA->left_end_screened  = false;
  WA->right_end_screened = false;

  WA->A_Olaps_For_Frag = 0;
  WA->B_Olaps_For_Frag = 0;

  Key = 0;
  for (j = 0;  j < G.Kmer_Len;  j ++)
    Key |= (uint64) (Bit_Equivalent [(int) * (P ++)]) << (2 * j);

  Sub = HASH_FUNCTION (Key);
  Shift = HASH_CHECK_FUNCTION (Key);
  Next_Key = (Key >> 2);
  Next_Key |= ((uint64) (Bit_Equivalent [(int) * P])) << (2 * (G.Kmer_Len - 1));
  Next_Sub = HASH_FUNCTION (Next_Key);
  Next_Shift = HASH_CHECK_FUNCTION (Next_Key);
  Next_Check = Hash_Check_Array [Next_Sub];

  if ((Hash_Check_Array [Sub] & (((Check_Vector_t) 1) << Shift)) != 0) {
    Ref = Hash_Find (Key, Sub, Window, & Where, & hi_hits);
    if (hi_hits) {
      WA->left_end_screened = true;
    }
    if (! getStringRefEmpty(Ref)) {
      while (true) {
        if (Frag_Num < getStringRefStringNum(Ref) + Hash_String_Num_Offset)
          Add_Ref  (Ref, Offset, WA);

        if (getStringRefLast(Ref))
          break;
        else {
          Ref = Extra_Ref_Space [++ Where];
          assert (! getStringRefEmpty(Ref));
        }
      }
    }
  }

  while ((* P) != '\0') {
    Window ++;
    Offset ++;

    Key = Next_Key;
    Shift = Next_Shift;
    Sub = Next_Sub;
    This_Check = Next_Check;
    P ++;
    Next_Key = (Key >> 2);
    Next_Key |= ((uint64)
                 (Bit_Equivalent [(int) * P])) << (2 * (G.Kmer_Len - 1));
    Next_Sub = HASH_FUNCTION (Next_Key);
    Next_Shift = HASH_CHECK_FUNCTION (Next_Key);
    Next_Check = Hash_Check_Array [Next_Sub];

    if ((This_Check & (((Check_Vector_t) 1) << Shift)) != 0) {
      Ref = Hash_Find (Key, Sub, Window, & Where, & hi_hits);
      if (hi_hits) {
        if (Offset < HOPELESS_MATCH) {
          WA->left_end_screened = true;
        }
        if (Frag_Len - Offset - G.Kmer_Len + 1 < HOPELESS_MATCH) {
          WA->right_end_screened = true;
        }
      }
      if (! getStringRefEmpty(Ref)) {
        while (true) {
          if (Frag_Num < getStringRefStringNum(Ref) + Hash_String_Num_Offset)
            Add_Ref  (Ref, Offset, WA);

          if (getStringRefLast(Ref))
            break;
          else {
            Ref = Extra_Ref_Space [++ Where];
            assert (! getStringRefEmpty(Ref));
          }
        }
      }
    }
  }


  Process_String_Olaps  (Frag, Frag_Len, Frag_Num, Dir, WA);
}

