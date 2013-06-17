
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

static const char *rcsid = "$Id: overlapInCore-Find_Overlaps.C,v 1.1 2011-07-30 01:16:05 brianwalenz Exp $";

#include "overlapInCore.H"

//  Add information for the match in  ref  to the list
//  starting at subscript  (* start) .  The matching window begins
//  offset  bytes from the beginning of this string.

static
void  Add_Match(String_Ref_t ref,
                int * start,
                int offset,
                int * consistent,
                Work_Area_t * wa) {
  int  * p, save;
  int  diag = 0, new_diag, expected_start = 0, num_checked = 0;
  int  move_to_front = FALSE;

  new_diag = getStringRefOffset(ref) - offset;

  for  (p = start;  (* p) != 0;  p = & (wa -> Match_Node_Space [(* p)] . Next))
    {
      expected_start = wa -> Match_Node_Space [(* p)] . Start
        + wa -> Match_Node_Space [(* p)] . Len
        - Kmer_Len + 1 + HASH_KMER_SKIP;
      //  Added  HASH_KMER_SKIP here --------------^^^^^^^^^^^^^^
      diag = wa -> Match_Node_Space [(* p)] . Offset
        - wa -> Match_Node_Space [(* p)] . Start;

      if  (expected_start < offset)
        break;
      if  (expected_start == offset)
        {
          if  (new_diag == diag)
            {
              wa -> Match_Node_Space [(* p)] . Len += 1 + HASH_KMER_SKIP;
              if  (move_to_front)
                {
                  save = (* p);
                  (* p) = wa -> Match_Node_Space [(* p)] . Next;
                  wa -> Match_Node_Space [save] . Next = (* start);
                  (* start) = save;
                }
              return;
            }
          else
            move_to_front = TRUE;
        }
      num_checked ++;
    }

  if  (wa -> Next_Avail_Match_Node == wa -> Match_Node_Size)
    {
      wa -> Match_Node_Size = (int32) (wa -> Match_Node_Size *
                                       MEMORY_EXPANSION_FACTOR);
      wa -> Match_Node_Space = (Match_Node_t *) safe_realloc
        (wa -> Match_Node_Space,
         wa -> Match_Node_Size * sizeof (Match_Node_t));
    }

  if  ((* start) != 0
       && (num_checked > 0
           || abs (diag - new_diag) > 3
           || offset < expected_start + Kmer_Len - 2))
    (* consistent) = FALSE;

  save = (* start);
  (* start) = wa -> Next_Avail_Match_Node;
  wa -> Next_Avail_Match_Node ++;

  wa -> Match_Node_Space [(* start)] . Offset = getStringRefOffset(ref);
  wa -> Match_Node_Space [(* start)] . Len = Kmer_Len;
  wa -> Match_Node_Space [(* start)] . Start = offset;
  wa -> Match_Node_Space [(* start)] . Next = save;

#if 0
  fprintf(stderr, "Add_Match()-- %3d offset %d len %d start %d next %d\n",
          *start,
          wa -> Match_Node_Space [(* start)] . Offset,
          wa -> Match_Node_Space [(* start)] . Len,
          wa -> Match_Node_Space [(* start)] . Start,
          wa -> Match_Node_Space [(* start)] . Next);
#endif
}



static void  Add_Ref
(String_Ref_t Ref, int Offset, Work_Area_t * WA)

//  Add information for  Ref  and all its matches to the global
//  hash table in  String_Olap_Space .  Grow the space if necessary
//  by  MEMORY_EXPANSION_FACTOR .  The matching window begins
//  Offset  bytes from the beginning of this string.

{
  uint32  Prev, StrNum, Sub;
  int  consistent;

  StrNum = getStringRefStringNum(Ref);
  Sub = (StrNum ^ (StrNum >> STRING_OLAP_SHIFT)) & STRING_OLAP_MASK;

  while  (WA -> String_Olap_Space [Sub] . Full
          && WA -> String_Olap_Space [Sub] . String_Num != StrNum)
    {
      Prev = Sub;
      Sub = WA -> String_Olap_Space [Sub] . Next;
      if  (Sub == 0)
        {
          if  (WA -> Next_Avail_String_Olap == WA -> String_Olap_Size)
            {
              WA -> String_Olap_Size = (int32) (WA -> String_Olap_Size *
                                                MEMORY_EXPANSION_FACTOR);
              WA -> String_Olap_Space = (String_Olap_t *) safe_realloc
                (WA -> String_Olap_Space,
                 WA -> String_Olap_Size * sizeof (String_Olap_t));
            }
          Sub = WA -> Next_Avail_String_Olap ++;
          WA -> String_Olap_Space [Prev] . Next = Sub;
          WA -> String_Olap_Space [Sub] . Full = FALSE;
          break;
        }
    }

  if  (! WA -> String_Olap_Space [Sub] . Full)
    {
      WA -> String_Olap_Space [Sub] . String_Num = StrNum;
      WA -> String_Olap_Space [Sub] . Match_List = 0;
      WA -> String_Olap_Space [Sub] . diag_sum = 0.0;
      WA -> String_Olap_Space [Sub] . diag_ct = 0;
      WA -> String_Olap_Space [Sub] . Next = 0;
      WA -> String_Olap_Space [Sub] . Full = TRUE;
      WA -> String_Olap_Space [Sub] . consistent = TRUE;
    }

  consistent = WA -> String_Olap_Space [Sub] . consistent;

  WA -> String_Olap_Space [Sub] . diag_sum += getStringRefOffset(Ref) - Offset;
  WA -> String_Olap_Space [Sub] . diag_ct ++;
  Add_Match (Ref, & (WA -> String_Olap_Space [Sub] . Match_List),
             Offset, & consistent, WA);

  WA -> String_Olap_Space [Sub] . consistent = consistent;

  return;
}




static String_Ref_t  Hash_Find
(uint64 Key, int64 Sub, char * S, int64 * Where, int * hi_hits)

//  Search for string  S  with hash key  Key  in the global
//  Hash_Table  starting at subscript  Sub .  Return the matching
//  reference in the hash table if there is one, or else a reference
//  with the  Empty bit set true.  Set  (* Where)  to the subscript in
//  Extra_Ref_Space  where the reference was found if it was found there.
//  Set  (* hi_hits)  to  TRUE  if hash table entry is found but is empty
//  because it was screened out, otherwise set to FALSE.

{
  String_Ref_t  H_Ref = 0;
  char  * T;
  unsigned char  Key_Check;
  int64  Ct, Probe;
  int  i;

  Key_Check = KEY_CHECK_FUNCTION (Key);
  Probe = PROBE_FUNCTION (Key);

  (* hi_hits) = FALSE;
  Ct = 0;
  do
    {
      for  (i = 0;  i < Hash_Table [Sub] . Entry_Ct;  i ++)
        if  (Hash_Table [Sub] . Check [i] == Key_Check)
          {
            int  is_empty;

            H_Ref = Hash_Table [Sub] . Entry [i];
            is_empty = getStringRefEmpty(H_Ref);
            if  (! getStringRefLast(H_Ref) && ! is_empty)
              {
                (* Where) = ((uint64)getStringRefStringNum(H_Ref) << OFFSET_BITS) + getStringRefOffset(H_Ref);
                H_Ref = Extra_Ref_Space [(* Where)];
              }
            T = Data + String_Start [getStringRefStringNum(H_Ref)] + getStringRefOffset(H_Ref);
            if  (strncmp (S, T, Kmer_Len) == 0)
              {
                if  (is_empty)
                  {
                    setStringRefEmpty(H_Ref, TRUELY_ONE);
                    (* hi_hits) = TRUE;
                  }
                return  H_Ref;
              }
          }
      if  (Hash_Table [Sub] . Entry_Ct < ENTRIES_PER_BUCKET)
        {
          setStringRefEmpty(H_Ref, TRUELY_ONE);
          return  H_Ref;
        }
      Sub = (Sub + Probe) % HASH_TABLE_SIZE;
    }  while  (++ Ct < HASH_TABLE_SIZE);

  setStringRefEmpty(H_Ref, TRUELY_ONE);
  return  H_Ref;
}






//  Find and output all overlaps and branch points between string
//   Frag  and any fragment currently in the global hash table.
//   Frag_Len  is the length of  Frag  and  Frag_Num  is its ID number.
//   Dir  is the orientation of  Frag .

void
Find_Overlaps (char Frag [], int Frag_Len, char quality [], AS_IID Frag_Num, Direction_t Dir, Work_Area_t * WA) {

  String_Ref_t  Ref;
  char  * P, * Window;
  uint64  Key, Next_Key;
  int64  Sub, Next_Sub, Where;
  Check_Vector_t  This_Check, Next_Check;
  int  Offset, Shift, Next_Shift;
  int  hi_hits;
  int  j;

  memset (WA -> String_Olap_Space, 0, STRING_OLAP_MODULUS * sizeof (String_Olap_t));
  WA -> Next_Avail_String_Olap = STRING_OLAP_MODULUS;
  WA -> Next_Avail_Match_Node = 1;

  assert (Frag_Len >= Kmer_Len);

  Offset = 0;
  P = Window = Frag;

  WA->left_end_screened  = FALSE;
  WA->right_end_screened = FALSE;

  WA -> A_Olaps_For_Frag = 0;
  WA -> B_Olaps_For_Frag = 0;

  Key = 0;
  for  (j = 0;  j < Kmer_Len;  j ++)
    Key |= (uint64) (Bit_Equivalent [(int) * (P ++)]) << (2 * j);

  Sub = HASH_FUNCTION (Key);
  Shift = HASH_CHECK_FUNCTION (Key);
  Next_Key = (Key >> 2);
  Next_Key |= ((uint64) (Bit_Equivalent [(int) * P])) << (2 * (Kmer_Len - 1));
  Next_Sub = HASH_FUNCTION (Next_Key);
  Next_Shift = HASH_CHECK_FUNCTION (Next_Key);
  Next_Check = Hash_Check_Array [Next_Sub];

  if  ((Hash_Check_Array [Sub] & (((Check_Vector_t) 1) << Shift)) != 0) {
    Ref = Hash_Find (Key, Sub, Window, & Where, & hi_hits);
    if  (hi_hits) {
      WA -> left_end_screened = TRUE;
    }
    if  (! getStringRefEmpty(Ref)) {
      while  (TRUE) {
        if  (Frag_Num < getStringRefStringNum(Ref) + Hash_String_Num_Offset)
          Add_Ref  (Ref, Offset, WA);

        if  (getStringRefLast(Ref))
          break;
        else {
          Ref = Extra_Ref_Space [++ Where];
          assert (! getStringRefEmpty(Ref));
        }
      }
    }
  }

  while  ((* P) != '\0') {
    Window ++;
    Offset ++;

    Key = Next_Key;
    Shift = Next_Shift;
    Sub = Next_Sub;
    This_Check = Next_Check;
    P ++;
    Next_Key = (Key >> 2);
    Next_Key |= ((uint64)
                 (Bit_Equivalent [(int) * P])) << (2 * (Kmer_Len - 1));
    Next_Sub = HASH_FUNCTION (Next_Key);
    Next_Shift = HASH_CHECK_FUNCTION (Next_Key);
    Next_Check = Hash_Check_Array [Next_Sub];

    if  ((This_Check & (((Check_Vector_t) 1) << Shift)) != 0) {
      Ref = Hash_Find (Key, Sub, Window, & Where, & hi_hits);
      if  (hi_hits) {
        if  (Offset < HOPELESS_MATCH) {
          WA -> left_end_screened = TRUE;
        }
        if  (Frag_Len - Offset - Kmer_Len + 1 < HOPELESS_MATCH) {
          WA -> right_end_screened = TRUE;
        }
      }
      if  (! getStringRefEmpty(Ref)) {
        while  (TRUE) {
          if  (Frag_Num < getStringRefStringNum(Ref) + Hash_String_Num_Offset)
            Add_Ref  (Ref, Offset, WA);

          if  (getStringRefLast(Ref))
            break;
          else {
            Ref = Extra_Ref_Space [++ Where];
            assert (! getStringRefEmpty(Ref));
          }
        }
      }
    }
  }


  Process_String_Olaps  (Frag, Frag_Len, quality, Frag_Num, Dir, WA);
}

