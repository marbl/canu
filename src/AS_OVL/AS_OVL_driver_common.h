
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
/*************************************************
* Module:  AS_OVL_driver.c
* Description:
*   Contains code to drive the overlapper, i.e., read messages, manage
*   stores and build threads
*************************************************/

/* RCS info
 * $Id: AS_OVL_driver_common.h,v 1.1.1.1 2004-04-14 13:52:27 catmandew Exp $
 * $Revision: 1.1.1.1 $
*/


#include  <unistd.h>

#include  "AS_OVL_delcher.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_fragStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_OVL_overlap.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_version.h"
/****vvvv*******  Screen Matches ***vvvv******/

VA_DEF(IntScreenMatch)

static VA_TYPE(IntScreenMatch) *ScreenMatches = NULL;

/****^^^^******  Screen Matches ****^^^^******/
static int64  First_Hash_Frag = -1;
static int64   Last_Hash_Frag;
static ReadStructp  myRead;
static int  Screen_Blocks_Used;
static int64  Total_Frags_Read = 0;
static int  Next_Distance_Index;
static int  Next_Fragment_Index;
static int  IID_Lo, IID_Hi;
static int  Frag_Segment_Lo;
static int  Frag_Segment_Hi;
#if  USE_THREADS
static pthread_mutex_t  Fragment_Range_Mutex;
#endif
static Batch_ID  Batch_Msg_UID = 0;
static IntBatch_ID  Batch_Msg_IID = 0;
static int  Batch_Num = 0;
static time_t  Now;


static void *  Choose_And_Process_Stream_Segment
    (void *);
static int  Choose_Hi_IID_Sub
    (uint32 List [], int lo, int n);
void  Cleanup_Work_Area
    (Work_Area_t * wa);
static void  Extract_Screen_Match_Info
    (ScreenedFragMesg * sfg_mesg, int numFragsRead);
static int  ReadFrags
    (int maxFrags, FragStoreHandle store, DistStore distStore, int argc, char **argv);


#if 0
#undef DEBUG
#define DEBUG 1
#undef MAX_HASH_STRINGS
#define MAX_HASH_STRINGS 20
#endif


static void  Check_VSize
    (void)
  {
   FILE  * fp;
   char  command [500];
   long unsigned  pid;
   static double  size, prev_size = 0.0;
   char  ch;

   pid = getpid ();
   sprintf (command, "ps -o vsize -p %lu > %lu.ps", pid, pid);
   system (command);

   sprintf (command, "%lu.ps", pid);
   fp = fopen (command, "r");
   assert (fp != NULL);

   fscanf (fp, "%s %lf%c", command, & size, & ch);
   if  (ch == 'K')
       size *= 1e-6;
   else if  (ch == 'M')
       size *= 1e-3;
   if  (size > prev_size)
       {
        fprintf (stderr, ">>> vsize increased to %.2fG from %.2fG\n",
                 size, prev_size);
        prev_size = size;
       }

   sprintf (command, "rm %lu.ps", pid);
   system (command);

   return;
  }



// **********************************************************************

int  OverlapDriver
    (int noOverlaps, int argc, char **argv)

//  This is the main control loop for the overlapper.
//  If  noOverlaps != 0 , then no overlaps are computed or output
//  (but the fragment store is created).

  {
   Frag_Stream  HashFragStream = 0;
   FragStore  NewFragStore;
   DistStore  NewDistStore;
#if  USE_THREADS
   pthread_attr_t  attr;
   pthread_t  * thread_id;
#endif
   Frag_Stream  * new_stream_segment;
   Frag_Stream  * old_stream_segment;
//   Work_Area_t  * driver_wa;
   Work_Area_t  * thread_wa;
   int64  first_new_frag = -1, last_new_frag = -1;
   int  i;

fprintf (stderr, "### sizeof (Work_Area_t) = " F_SIZE_T "\n",
         sizeof (Work_Area_t));
#if  USE_THREADS
fprintf (stderr, "### Using %d pthreads  %d hash bits  %d bucket entries\n",
         Num_PThreads, HASH_MASK_BITS, ENTRIES_PER_BUCKET);
#else
fprintf (stderr, "### Not using threads  %d hash bits  %d bucket entries\n",
         HASH_MASK_BITS, ENTRIES_PER_BUCKET);
#endif

#if  USE_THREADS
   thread_id = (pthread_t *) Safe_calloc
                   (Num_PThreads, sizeof (pthread_t));
#endif

   new_stream_segment = (Frag_Stream *) Safe_calloc
                   (Num_PThreads, sizeof (Frag_Stream));
   old_stream_segment = (Frag_Stream *) Safe_calloc
                   (Num_PThreads, sizeof (Frag_Stream));
//   driver_wa = (Work_Area_t *) Safe_malloc (sizeof (Work_Area_t));
   thread_wa = (Work_Area_t *) Safe_calloc
                   (Num_PThreads, sizeof (Work_Area_t));
   if  (Contig_Mode)
       NewFragStore
           = createFragStore (NULL, "temp", 1);
   else if  (! LSF_Mode)
       {
        first_new_frag = getLastElemFragStore (OldFragStore) + 1;
        NewFragStore
            = createFragStore (NULL, "temp", first_new_frag);
       }
   if  (! LSF_Mode)
       NewDistStore
           = createDistStore (NULL, getLastElemStore (OldDistStore) + 1);
   if  (! LSF_Mode)
       HashFragStream = openFragStream (NewFragStore, NULL, 0);
   for  (i = 0;  i < Num_PThreads;  i ++)
     {
      if  (! LSF_Mode)
          new_stream_segment [i] = openFragStream (NewFragStore, NULL, 0);
      old_stream_segment [i] = openFragStream (OldFragStore, NULL, 0);
     }

   if  (noOverlaps == 0)
       {
#if  USE_THREADS
        if  (Num_PThreads > 1)
            {
             pthread_attr_init (& attr);
             pthread_attr_setstacksize (& attr, THREAD_STACKSIZE);
             pthread_mutex_init (& Fragment_Range_Mutex, NULL);
             pthread_mutex_init (& FragStore_Mutex, NULL);
             pthread_mutex_init (& Write_Proto_Mutex, NULL);
             pthread_mutex_init (& Log_Msg_Mutex, NULL);
            }
#endif
//        Initialize_Work_Area (driver_wa, -1);
        Initialize_Work_Area (thread_wa, 0);
        for  (i = 1;  i < Num_PThreads;  i ++)
          Initialize_Work_Area (thread_wa + i, i);
       }

   myRead = new_ReadStruct ();

#if  SHOW_PROGRESS
Start_Time = clock ();
#endif

   if  (Contig_Mode)
       Next_Fragment_Index = 1;
     else
       Next_Fragment_Index = getLastElemFragStore (OldFragStore) + 1;
   if  (! LSF_Mode)
       Next_Distance_Index = getLastElemStore (OldDistStore) + 1;


//  Call  readFrags  until we've exhausted the input stream.
//  readFrags  reads from input, and adds the fragments to
//  NewFragStore .

#if  USE_SOURCE_FIELD
Source_Log_File = File_Open ("ovl-srcinfo.log", "w");
#endif

   if  (LSF_Mode)
       {
        int  id;

        if  (Contig_Mode)
            id = getFirstElemFragStore (BACtigStore);
          else
            id = getFirstElemFragStore (OldFragStore);
        if  (Lo_Hash_Frag < id)
            Lo_Hash_Frag = id;

        if  (Contig_Mode)
            id = getLastElemFragStore (BACtigStore);
          else
            id = getLastElemFragStore (OldFragStore);
        if  (id < Hi_Hash_Frag)
            Hi_Hash_Frag = id;
       }

   while (ReadFrags (MAX_HASH_STRINGS, NewFragStore, NewDistStore, argc, argv))
     {
      int startIndex;

      if  (noOverlaps == 0)
          {
           FragStore  curr_frag_store;
           FragStore  hash_frag_store;
           int  highest_old_frag, lowest_old_frag;
           int  status;

           if  (LSF_Mode)
               {
                if  (Contig_Mode)
                    {
                     hash_frag_store
                         = loadFragStorePartial (BACtig_Store_Path,
                                                 First_Hash_Frag,
                                                 Last_Hash_Frag);
                     fprintf (stderr,
                              "loadFragStorePartial  first = " F_S64 "  last = " F_S64 "\n",
                              First_Hash_Frag, Last_Hash_Frag);
                     assert (0 < First_Hash_Frag
                               && First_Hash_Frag <= Last_Hash_Frag
                               && Last_Hash_Frag
                                    <= getLastElemFragStore (BACtigStore));
                    }
                  else
                    {
                     hash_frag_store
                         = loadFragStorePartial (Frag_Store_Path,
                                                 First_Hash_Frag,
                                                 Last_Hash_Frag);
                     fprintf (stderr,
                              "loadFragStorePartial  first = " F_S64 "  last = " F_S64 "\n",
                              First_Hash_Frag, Last_Hash_Frag);
                     assert (0 < First_Hash_Frag
                               && First_Hash_Frag <= Last_Hash_Frag
                               && Last_Hash_Frag
                                    <= getLastElemFragStore (OldFragStore));
                    }
                HashFragStream = openFragStream (hash_frag_store, NULL, 0);
                resetFragStream (HashFragStream, First_Hash_Frag,
                                 Last_Hash_Frag);
                startIndex = First_Hash_Frag;
               }
             else
               {
                resetFragStream (HashFragStream, STREAM_FROMSTART,
                                 STREAM_UNTILEND);
                startIndex = getFirstElemFragStore (NewFragStore);
               }


/* Create the hash table from the HashFragStream */

           Build_Hash_Index (HashFragStream, startIndex, myRead);

           if  (Last_Hash_Frag_Read < Last_Hash_Frag)
               {
                fprintf (stderr, "!!! Hash table did not read all frags\n");
                fprintf (stderr, "    Read " F_U32 " instead of " F_S64 "\n",
                         Last_Hash_Frag_Read, Last_Hash_Frag);
                Last_Hash_Frag = Last_Hash_Frag_Read;
               }

#if  DO_KMER_HITS_PROFILE
break;
#endif

           resetFragStream (HashFragStream, STREAM_FROMSTART,
                            STREAM_UNTILEND);

           if  (! LSF_Mode && ! Contig_Mode)
               {
                /* Stream HashFragStream by the hash table, to compute overlaps */

                Frag_Segment_Lo = getFirstElemFragStore (NewFragStore);
                Frag_Segment_Hi = getLastElemFragStore (NewFragStore);
                curr_frag_store = NewFragStore;

#if  USE_THREADS
                for  (i = 1;  i < Num_PThreads;  i ++)
                  {
                   thread_wa [i] . stream_segment = new_stream_segment [i];
                   status = pthread_create
                                (thread_id + i, & attr,
                                 Choose_And_Process_Stream_Segment,
                                 thread_wa + i);
                   if  (status != 0)
                       {
                        fprintf (stderr, "pthread_create error at line %d:  %s\n",
                                 __LINE__, strerror (status));
                        exit (-3);
                       }
                  }
#endif

                thread_wa [0] . stream_segment = new_stream_segment [0];
                Choose_And_Process_Stream_Segment (thread_wa);

#if  USE_THREADS
                for  (i = 1;  i < Num_PThreads;  i ++)
                  {
                   void  * ptr;

                   status = pthread_join  (thread_id [i], & ptr);
                   if  (status != 0)
                       {
                        fprintf (stderr, "pthread_join error at line %d:  %s\n",
                                 __LINE__, strerror (status));
                        exit (-3);
                       }
                  }
#endif
               }


           lowest_old_frag = getFirstElemFragStore (OldFragStore);
           highest_old_frag = getLastElemFragStore (OldFragStore);
           if  (lowest_old_frag < Lo_Old_Frag)
               lowest_old_frag = Lo_Old_Frag;
           if  (highest_old_frag > Hi_Old_Frag)
               highest_old_frag = Hi_Old_Frag;
           if  (LSF_Mode && ! Contig_Mode
                  && highest_old_frag > Last_Hash_Frag)
               highest_old_frag = Last_Hash_Frag;
           if  (IID_List != NULL)
               {
                if  (lowest_old_frag < IID_List [0])
                    lowest_old_frag = IID_List [0];
                if  (highest_old_frag > IID_List [IID_List_Len - 1])
                    highest_old_frag = IID_List [IID_List_Len - 1];
                IID_Lo = 0;
               }

#if  ! SCREEN_CHECK_ONLY
           while  (lowest_old_frag <= highest_old_frag)
             {
              Frag_Segment_Lo = lowest_old_frag;
              if  (IID_List == NULL)
                  {
                   Frag_Segment_Hi = Frag_Segment_Lo + MAX_FRAGS_IN_MEMORY_STORE - 1;
                   if  (Frag_Segment_Hi > highest_old_frag)
                       Frag_Segment_Hi = highest_old_frag;
                  }
                else
                  {
                   IID_Hi = Choose_Hi_IID_Sub (IID_List, IID_Lo, IID_List_Len);
                   Frag_Segment_Hi = IID_List [IID_Hi];
                  }

              curr_frag_store
                  = loadFragStorePartial
                        (Frag_Store_Path, Frag_Segment_Lo, Frag_Segment_Hi);
              fprintf (stderr,
                       "loadFragStorePartial  first = %d  last = %d\n",
                       Frag_Segment_Lo, Frag_Segment_Hi);
              assert (0 < Frag_Segment_Lo
                        && Frag_Segment_Lo <= Frag_Segment_Hi
                        && Frag_Segment_Hi
                             <= getLastElemFragStore (OldFragStore));

              for  (i = 0;  i < Num_PThreads;  i ++)
                {
                 old_stream_segment [i] = openFragStream (curr_frag_store, NULL, 0);
                 resetFragStream (old_stream_segment [i], Frag_Segment_Lo,
                                  Frag_Segment_Hi);
                }


              Now = time (NULL);
              fprintf (stderr, "### starting old fragments   %s\n", ctime (& Now));
#if  USE_THREADS
              for  (i = 1;  i < Num_PThreads;  i ++)
                {
                 thread_wa [i] . stream_segment = old_stream_segment [i];
                 status = pthread_create
                              (thread_id + i, & attr,
                               Choose_And_Process_Stream_Segment,
                               thread_wa + i);
                 if  (status != 0)
                     {
                      fprintf (stderr, "pthread_create error at line %d:  %s\n",
                               __LINE__, strerror (status));
                      exit (-3);
                     }
                }
#endif

              thread_wa [0] . stream_segment = old_stream_segment [0];
              Choose_And_Process_Stream_Segment (thread_wa);

#if  USE_THREADS
              for  (i = 1;  i < Num_PThreads;  i ++)
                {
                 void  * ptr;

                 status = pthread_join  (thread_id [i], & ptr);
                 if  (status != 0)
                     {
                      fprintf (stderr, "pthread_join error at line %d:  %s\n",
                               __LINE__, strerror (status));
                      exit (-3);
                     }
                }
#endif

              Now = time (NULL);
              fprintf (stderr, "### done old fragments   %s", ctime (& Now));
              for  (i = 0;  i < Num_PThreads;  i ++)
                closeFragStream (old_stream_segment [i]);
              closeFragStore (curr_frag_store);

              if  (IID_List == NULL)
                  lowest_old_frag += MAX_FRAGS_IN_MEMORY_STORE;
                else
                  {
                   IID_Lo = IID_Hi + 1;
                   if  (IID_Lo < IID_List_Len)
                       lowest_old_frag = IID_List [IID_Lo];
                     else
                       lowest_old_frag = INT_MAX;
                  }
             }
#endif
                 
           if  (LSF_Mode)
               {
                closeFragStream (HashFragStream);
                closeFragStore (hash_frag_store);
               }
          }

      if  (! LSF_Mode)
          {
           if  (Contig_Mode)
               last_new_frag = getLastElemFragStore (NewFragStore);
             else
               {
                concatFragStore (OldFragStore, NewFragStore);
                commitFragStore (OldFragStore);
                last_new_frag = getLastElemFragStore (OldFragStore);

                /* If we've added dst records, save them persistently */

                if  (getLastElemStore (NewDistStore) > getLastElemStore (OldDistStore))
                    {
                     concatStore (OldDistStore, NewDistStore);
                    }
                commitStore (OldDistStore);
                resetDistStore (NewDistStore, getLastElemStore (OldDistStore) + 1);
               }

           resetFragStore (NewFragStore, last_new_frag + 1);
          }

#if  SHOW_PROGRESS
Stop_Time = clock ();
fprintf (stderr, "Table:%d %7.1f sec %7ld olaps\n",
         Table_Ct, (double) (Stop_Time - Start_Time) / CLOCKS_PER_SEC,
         Olap_Ct);
Table_Ct ++;
Olap_Ct = 0;
Start_Time = clock ();
#endif

      Now = time (NULL);
      fprintf (stderr, "### Done batch #%d   %s", Batch_Num, ctime (& Now));

#if  ANALYZE_HITS && ! DO_KMER_HITS_PROFILE
Output_High_Hit_Frags ();
#endif
     }

#if  DO_KMER_HITS_PROFILE
Profile_Hits ();
#endif

   if  (first_new_frag >= 0)
       {
        FILE  * fp = File_Open ("new.range", "w");

        fprintf (fp, "batch:" F_UID "\n", Batch_Msg_UID);
        fprintf (fp, "lo:" F_S64 "  hi:" F_S64 "\n", first_new_frag, last_new_frag);
        fclose (fp);
       }

   /* Handle the pathological case where we read ONLY distance records */
   /* If we've added dst records, save them persistently */

   if  (! LSF_Mode)
       {
        if  (getLastElemStore (NewDistStore) > getLastElemStore (OldDistStore))
          {
            concatStore (OldDistStore, NewDistStore);
          }
        commitStore (OldDistStore);
       }

//   Cleanup_Work_Area (driver_wa);
   Cleanup_Work_Area (thread_wa);
   for  (i = 1;  i < Num_PThreads;  i ++)
     Cleanup_Work_Area (thread_wa + i);
//   free (driver_wa);
   free (thread_wa);
#if  USE_THREADS
   free (thread_id);
#endif
   free (new_stream_segment);
   free (old_stream_segment);

#if  ANALYZE_HITS
fclose (High_Hits_File);
#endif

#if  USE_SOURCE_FIELD
fclose (Source_Log_File);
#endif

   closeFragStream (HashFragStream);
   if  (! LSF_Mode)
       {
        closeDistStore (OldDistStore);
        closeDistStore (NewDistStore);
        closeFragStore (NewFragStore);
       }

   fprintf (stderr, "Total fragments read = " F_S64 "\n", Total_Frags_Read);

   return  0;
  }



/******************************************************************************/

static int  Choose_Hi_IID_Sub
    (uint32 List [], int lo, int n)

//  Return subscript in range  lo .. n  so that fragment IID's in
//  that range of  List  are suitable for processing in a round of
//  overlaps.

  {
   int  i;

   for  (i = lo + 1;
           i < n
             && List [i] - List [i - 1] < IID_GAP_LIMIT
             && List [i] - List [lo] < MAX_FRAGS_IN_MEMORY_STORE;
           i ++)
     ;

   return  i - 1;
  }



/******************************************************************************/

void  Cleanup_Work_Area
    (Work_Area_t * wa)

//  Free memory allocated in  (* wa) .

  {
   free (wa -> String_Olap_Space);
   free (wa -> Match_Node_Space);
   delete_ReadStruct (wa -> myRead);

   return;
  }



/******************************************************************************/

void  Coalesce_Screen_Info
    (Screen_Range_t space [], int lo, int * top)

//  Sort entries in  space [lo .. (* top) - 1]
//  Then combine overlapping adjacent entries.  Reduce  (* top)
//  by the number of combined entries.

  {
   int  i, j, hi;

   hi = (* top - 1);

   if  (hi <= lo)
       return;

   //  Sort ascending by  bgn  value
   for  (i = lo;  i < hi;  i ++)
     for  (j = i + 1;  j <= hi;  j ++)
       if  (space [i] . bgn > space [j] . bgn)
           {
            Screen_Range_t  save = space [i];

            space [i] = space [j];
            space [j] = save;
           }

   //  Combine overlapping entries.  Regard as overlapping if gap between
   //  is too short to allow a k-mer hit
   for  (i = lo, j = i + 1;  j <= hi;  j ++)
     {
      if  (space [j] . bgn >= space [i] . end + WINDOW_SIZE - WINDOW_SCREEN_OLAP)
          space [++ i] = space [j];
        else if  (space [j] . end
                    > space [i] . end)
          space [i] . end = space [j] . end;
     }

   (* top) = i + 1;

   return;
  }



/******************************************************************************/

static void *  Choose_And_Process_Stream_Segment
    (void * ptr)

//  Find all overlaps between frags  Frag_Segment_Lo .. Frag_Segment_Hi
//  in the stream in  ptr  and the frags in the hash table.

  {
   Work_Area_t  * WA = (Work_Area_t *) (ptr);

   while  (TRUE)
     {
      int  lo, hi;

#if  USE_THREADS
      if  (Num_PThreads > 1)
          pthread_mutex_lock (& Fragment_Range_Mutex);
#endif

      if  (IID_List == NULL)
          {
           lo = Frag_Segment_Lo;
           Frag_Segment_Lo += MAX_FRAGS_PER_THREAD;
           hi = Frag_Segment_Lo - 1;
          }
        else
          {
           if  (IID_Lo > IID_Hi)
               lo = hi = INT_MAX;
             else
               lo = hi = IID_List [IID_Lo ++];
          }

#if  USE_THREADS
      if  (Num_PThreads > 1)
          pthread_mutex_unlock (& Fragment_Range_Mutex);
#endif

      if  (IID_List != NULL && hi > Frag_Segment_Hi)
          break;
      if  (hi > Frag_Segment_Hi)
          hi = Frag_Segment_Hi;
      if  (lo > hi)
          break;
      
#if  USE_THREADS
      if  (Num_PThreads > 1)
          pthread_mutex_lock (& FragStore_Mutex);
#endif

      resetFragStream (WA -> stream_segment, lo, hi);

#if  USE_THREADS
      if  (Num_PThreads > 1)
          pthread_mutex_unlock (& FragStore_Mutex);
#endif

      Process_Overlaps (WA -> stream_segment, WA);
     }

#if  USE_THREADS
   if  (Num_PThreads > 1 && WA -> thread_id > 0)
       pthread_exit (ptr);
#endif

   return  ptr;
  }



/*************************************************************************/

static void  Extract_Screen_Match_Info
    (ScreenedFragMesg * sfg_mesg, int numFragsRead)

//  Extract the screened regions specified in  (* sfg_mesg)  and
//  save them (with any overlaps combined) in the global  Screen_Space .
//  Make  Screen_Sub [numFragsRead]  have the offset to the entries.

  {
   IntScreenMatch  * p;
   int  match_ct;
   int  lo = Screen_Blocks_Used;

   if(ScreenMatches == NULL){
     ScreenMatches = CreateVA_IntScreenMatch(512);
   }else{
     ResetIntScreenMatch(ScreenMatches);
   }

   for  (p = sfg_mesg -> screened;  p != NULL;  p = p -> next)
     {
     AppendIntScreenMatch(ScreenMatches, p);
     if  (p -> relevance & AS_OVL_HEED_RPT)
         {
          SeqInterval  s;

          s = p -> where;
          if  (Screen_Blocks_Used >= Screen_Space_Size)
              {
               Screen_Space_Size *= MEMORY_EXPANSION_FACTOR;
               fprintf (stderr,
                        "### realloc  Screen_Space  Screen_Space_Size = %d\n",
                        Screen_Space_Size);
               Screen_Space = (Screen_Range_t *) Safe_realloc
                                  (Screen_Space,
                                  Screen_Space_Size * sizeof (Screen_Range_t));
              }
          if  (s . bgn < 0)
              {
               fprintf
                   (stderr,
                    "ERROR:  Screen begin = %d on fragment %d;  reset to 0\n",
                    s . bgn, sfg_mesg -> iaccession);
               s . bgn = 0;
              }
          Screen_Space [Screen_Blocks_Used] . bgn = s . bgn;
          if  (s . end > AS_READ_MAX_LEN)
              {
               fprintf
                   (stderr,
                    "ERROR:  Screen end = %d on fragment %d;  reset\n",
                    s . end, sfg_mesg -> iaccession);
               s . end = AS_READ_MAX_LEN;
              }
          Screen_Space [Screen_Blocks_Used] . end = s . end;
          Screen_Space [Screen_Blocks_Used] . last = FALSE;
          Screen_Blocks_Used ++;
         }
     }

   // Make  next  pointers point to entries in  ScreenMatches
   // Can't do earlier since  ScreenMatches  might be realloc'd.

   match_ct = GetNumIntScreenMatchs (ScreenMatches);
   if  (match_ct > 0)
       {
        IntScreenMatch  * p = GetIntScreenMatch (ScreenMatches, 0);
        IntScreenMatch  * q;
        int  i;
        
        for  (i = 1;  i < match_ct;  i ++)
          {
           q = GetIntScreenMatch (ScreenMatches, i);
           p -> next = q;
           p = q;
          }

        p -> next = NULL;
       }

   if  (lo == Screen_Blocks_Used)
       Screen_Sub [numFragsRead - 1] = 0;
     else
       {
        Coalesce_Screen_Info (Screen_Space, lo, & Screen_Blocks_Used);
        Screen_Sub [numFragsRead - 1] = lo;
        Screen_Space [Screen_Blocks_Used - 1] . last = TRUE;
       }

   return;
  }



/******************************************************************************/

/* function ReadFrags:
   Read up to maxFrags, but maybe a little less.

 */
static int  ReadFrags
    (int maxFrags, FragStoreHandle store, DistStore distStore, int argc, char **argv)

  {
  int numFragsRead = 0;
  MessageType imesgtype;
  GenericMesg   *pmesg;
  InternalDistMesg  *idt_mesg;
  ScreenedFragMesg *sfg_mesg, hack_sfg_mesg;
  InternalFragMesg *ifg_mesg;
  OFGMesg ofg_mesg;
  int  first_distance, first_hash_frag;
  int64  total_len = 0;

  Now = time (NULL);
  fprintf (stderr, "### Start batch #%d   %s", ++ Batch_Num, ctime (& Now));

  // In LSF mode there is no input file.  Instead, the hash table is built
  // from fragments from the fragment store.  In this case, we just
  // set values of global variables  First_Hash_Frag  and  Last_Hash_Frag
  // to indicate the range of fragments to get from that store.

  if  (LSF_Mode)
      {
       if  (First_Hash_Frag == -1)    // First time through
           First_Hash_Frag = Lo_Hash_Frag;
         else
           First_Hash_Frag = Last_Hash_Frag + 1;

       if  (First_Hash_Frag > Hi_Hash_Frag)
           return  FALSE;

       {
	 int64 temp = First_Hash_Frag + maxFrags - 1;
	 Last_Hash_Frag = ( temp < Hi_Hash_Frag ? temp : Hi_Hash_Frag );
       }
       return  TRUE;
      }


  Screen_Blocks_Used = 1;
  first_distance = Next_Distance_Index;
  first_hash_frag = Next_Fragment_Index;

  /* Read maxFrags fragments from the stream */
  while  ((numFragsRead < maxFrags)
            && HASH_EXPANSION_FACTOR * total_len
                   <= HASH_TABLE_SIZE * ENTRIES_PER_BUCKET
            && (EOF != Read_Msg_Fn (In_Stream, & pmesg)))
    {
    imesgtype = pmesg->t;
    switch(imesgtype){
    case MESG_IDT:
      {
	DistRecord distRecord;
	idt_mesg = (InternalDistMesg*) pmesg->m;
#if DEBUG
	fprintf(stderr,"Read DST message %c ( " F_S64 ", %d ) Next_Distance_Index = %d\n", 
		idt_mesg->action,
		idt_mesg->eaccession,
		idt_mesg->iaccession,
		Next_Distance_Index);
#endif
	distRecord . UID = idt_mesg -> eaccession;
	distRecord . IID = idt_mesg -> iaccession;
        distRecord . mean = idt_mesg -> mean;
        distRecord . stddev = idt_mesg -> stddev;
        distRecord . deleted = 0;
        distRecord . spare = 0;

	switch(idt_mesg->action){
	case AS_ADD:
	  {
	    if(distRecord.IID != Next_Distance_Index){
	      fprintf(stderr,"*** Fatal Error -- distance record should have IID %d not %d\n",
		      Next_Distance_Index, distRecord.IID);
	      exit(1);
	    }

	    Next_Distance_Index++;

	    //	    fprintf(stderr,"* Distance Message (" F_S64 ",%d)\n", 
	    //    distRecord.UID, distRecord.IID);

	    appendDistStore(distStore, &distRecord);
	  }
	  break;
	case AS_DELETE:
	  if  (distRecord . IID > Next_Distance_Index)
              {
	       fprintf (stderr, "*** Fatal Error -- distance record %d doesn't exist\n",
		        idt_mesg -> iaccession);
	       exit (EXIT_FAILURE);
              }
          else if  (idt_mesg -> iaccession < first_distance)
              deleteDistStore (OldDistStore, idt_mesg -> iaccession);
            else
              deleteDistStore (distStore, idt_mesg -> iaccession);
          
          break;

        case  AS_REDEFINE :
	  if  (distRecord . IID >= Next_Distance_Index)
              {
	       fprintf (stderr,
                        "*** Fatal Error -- Attempt to redefine distance record %d\n"
                        "***   which has not yet been seen\n",
		        idt_mesg -> iaccession);
	       exit (EXIT_FAILURE);
              }
          else if  (idt_mesg -> iaccession < first_distance)
              setDistStore (OldDistStore, idt_mesg -> iaccession, & distRecord);
            else
              setDistStore (distStore, idt_mesg -> iaccession, & distRecord);

          break;
        default:
          assert(0);
        }

#if  ! (FOR_CARL_FOSLER || SHOW_SNPS)
       Write_Msg_Fn (Out_Stream,pmesg);
#endif
      }
      break;

    case MESG_IFG:
      ifg_mesg = (InternalFragMesg*) pmesg->m;
      Transfer_IFG_to_SFG_AS (ifg_mesg, &hack_sfg_mesg);
      pmesg->m = &hack_sfg_mesg;
      pmesg->t = MESG_SFG;
      /*** FALL THROUGH ***/
    case MESG_SFG:
        {
          /* Put the record where it belongs in the array.
             This array is indexed by the overlaps. */
          int  clear_len;

          sfg_mesg = (ScreenedFragMesg*) pmesg->m;
          Transfer_SFG_to_OFG_AS (sfg_mesg, &ofg_mesg);
          pmesg->m = &ofg_mesg;
          pmesg->t = MESG_OFG;

#if DEBUG
        fprintf(stderr,"Read IFG/SFG message %c ( " F_S64 ", %d ) \n",
                idt_mesg->action,
                sfg_mesg->eaccession,
                sfg_mesg->iaccession);
#endif
          switch(sfg_mesg->action){
          case AS_ADD:
            numFragsRead++;
            clear_len = ofg_mesg.clear_rng.end - ofg_mesg.clear_rng.bgn;
#if DEBUG
            fprintf(stderr,"Read message %c (" F_S64 ", %d) %d\n", 
                    sfg_mesg->action, sfg_mesg->eaccession,
                    sfg_mesg->iaccession, Next_Fragment_Index);
#endif
            if(sfg_mesg->iaccession != Next_Fragment_Index){
              fprintf(stderr,"*** Fatal Error -- fragment record should have IID %d not %d\n",
                      Next_Fragment_Index, sfg_mesg->iaccession);
              exit(1);
            }


            if  (sfg_mesg -> screened != NULL)
                Extract_Screen_Match_Info
                    (sfg_mesg, numFragsRead);
              else
                Screen_Sub [numFragsRead - 1] = 0;
                

            Next_Fragment_Index++;

            /* Add it to the fragment Store, as well */
            setAccID_ReadStruct(myRead, ofg_mesg.eaccession);
            setReadIndex_ReadStruct(myRead, ofg_mesg.iaccession);
            setReadType_ReadStruct(myRead, ofg_mesg.type);
            stripWhiteSpace(Sequence_Buffer, sfg_mesg->sequence, AS_READ_MAX_LEN * 2);
            stripWhiteSpace(Quality_Buffer, sfg_mesg->quality, AS_READ_MAX_LEN * 2);
            setSequence_ReadStruct(myRead, Sequence_Buffer, Quality_Buffer);
            setSource_ReadStruct(myRead, sfg_mesg->source);
            setEntryTime_ReadStruct(myRead, ofg_mesg.entry_time);
	    // These clear ranges are the original ones, not OVL-modified.
            setClearRegion_ReadStruct
                (myRead, ofg_mesg.clear_rng.bgn, ofg_mesg.clear_rng.end,
		 READSTRUCT_ORIGINAL);  
            setLocalePos_ReadStruct
                (myRead,ofg_mesg.locale_pos.bgn, ofg_mesg.locale_pos.end);
            // changed by Knut Reinert
            // due to gatekeeper changes
            setLocID_ReadStruct(myRead,ofg_mesg.ilocale);
            if  (sfg_mesg -> screened != NULL)
            {
              IntScreenMatch *matches = GetIntScreenMatch(ScreenMatches,0);
              setScreenMatches_ReadStruct(myRead, GetNumIntScreenMatchs(ScreenMatches), matches);
            }else{
              setScreenMatches_ReadStruct(myRead, 0, NULL);
            }
            total_len += clear_len;

            appendFragStore(store, myRead);

          break;

        case AS_DELETE:
          if  (sfg_mesg -> iaccession < first_hash_frag)
              deleteFragStore (OldFragStore, sfg_mesg -> iaccession);
            else
              deleteFragStore (store, sfg_mesg -> iaccession);
	  break;

	default:
	  assert(0);
	  break;
	  }

#if  ! (FOR_CARL_FOSLER || SHOW_SNPS)
	Write_Msg_Fn (Out_Stream, pmesg);
#endif
	}
	break;

      case MESG_ADT:
	{
	  AuditMesg *adt_mesg;

	  adt_mesg = (AuditMesg*) pmesg->m;
	  pmesg->t = MESG_ADT;

	  VersionStampADT(adt_mesg, argc, argv);

#if  ! (FOR_CARL_FOSLER || SHOW_SNPS)
	Write_Msg_Fn (Out_Stream, pmesg);
#endif
	}
	break;

      case  MESG_IBA :
        {
         InternalBatchMesg  * iba_mesg;

         iba_mesg = (InternalBatchMesg*) pmesg -> m;
         Batch_Msg_UID = iba_mesg -> eaccession;
         Batch_Msg_IID = iba_mesg -> iaccession;

         // fall through
        }

      case MESG_ILK:
      case MESG_ISN:
      case MESG_RPT:
      case MESG_IRP :
      case MESG_IBC:

#if DEBUG
	fprintf(stderr,"Read ILK/IJN/ISN/RPT\n");
#endif

#if  ! (FOR_CARL_FOSLER || SHOW_SNPS)
	Write_Msg_Fn (Out_Stream, pmesg);
#endif
	break;

      default:
	fprintf(stderr,"* Oops: Read Message with type imesgtype = %d\n", imesgtype);
	WriteProtoMesg_AS(stderr,pmesg);      
	
	exit(1);
	}
      }

  fprintf (stderr, "### %d fragments read\n", numFragsRead);

  Total_Frags_Read += numFragsRead;

  return(numFragsRead);
}


/******************************************************************************/

/* Function stripWhiteSpace:
   Input:  source   string of maximum length maxlen
           maxlen   maximum length of source
   Output: target
   
   Description:
     Copy non-space characters from source to target.
*/

void  stripWhiteSpace
    (char *target, char *source, int maxlen)

  {
  int i = 0;
  *target = '\0';
  while(i < maxlen){
    if(!isspace(*source)){
      *target++ = *source;
      i++;
    }
    if(*source == '\0')
      break;
    source++;
  }

}
