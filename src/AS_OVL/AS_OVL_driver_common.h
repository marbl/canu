
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
 * $Id: AS_OVL_driver_common.h,v 1.16 2007-02-22 14:44:40 brianwalenz Exp $
 * $Revision: 1.16 $
*/


#include  <unistd.h>

#include  "AS_OVL_delcher.h"
#include  "AS_PER_gkpStore.h"
#include  "AS_PER_genericStore.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_OVL_overlap.h"
#include  "AS_UTL_Var.h"
#include  "AS_UTL_version.h"

static int64  First_Hash_Frag = -1;
static int64   Last_Hash_Frag;
static fragRecord  *myRead;
static int  Screen_Blocks_Used;
static int64  Total_Frags_Read = 0;
static int  Next_Distance_Index;
static int  Next_Fragment_Index;
static int  IID_Lo, IID_Hi;
static int  Frag_Segment_Lo;
static int  Frag_Segment_Hi;
static pthread_mutex_t  Fragment_Range_Mutex;
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
static int  ReadFrags
    (int maxFrags);


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
   FragStream  *HashFragStream = NULL;
   pthread_attr_t  attr;
   pthread_t  * thread_id;
   FragStream **new_stream_segment;
   FragStream **old_stream_segment;
//   Work_Area_t  * driver_wa;
   Work_Area_t  * thread_wa;
   int64  first_new_frag = -1, last_new_frag = -1;
   int  i;

fprintf (stderr, "### sizeof (Work_Area_t) = " F_SIZE_T "\n",
         sizeof (Work_Area_t));
fprintf (stderr, "### Using %d pthreads  %d hash bits  %d bucket entries\n",
         Num_PThreads, Hash_Mask_Bits, ENTRIES_PER_BUCKET);

   thread_id = (pthread_t *) safe_calloc
                   (Num_PThreads, sizeof (pthread_t));

   new_stream_segment = (FragStream **) safe_calloc
                   (Num_PThreads, sizeof (FragStream *));
   old_stream_segment = (FragStream **) safe_calloc
                   (Num_PThreads, sizeof (FragStream *));
//   driver_wa = (Work_Area_t *) safe_malloc (sizeof (Work_Area_t));
   thread_wa = (Work_Area_t *) safe_calloc
                   (Num_PThreads, sizeof (Work_Area_t));

   for  (i = 0;  i < Num_PThreads;  i ++)
     {
      old_stream_segment [i] = openFragStream (OldFragStore, FRAG_S_INF | FRAG_S_SEQ | FRAG_S_QLT);
     }

   if  (noOverlaps == 0)
       {
        if  (Num_PThreads > 1)
            {
             pthread_attr_init (& attr);
             pthread_attr_setstacksize (& attr, THREAD_STACKSIZE);
             pthread_mutex_init (& Fragment_Range_Mutex, NULL);
             pthread_mutex_init (& FragStore_Mutex, NULL);
             pthread_mutex_init (& Write_Proto_Mutex, NULL);
             pthread_mutex_init (& Log_Msg_Mutex, NULL);
            }
        Initialize_Work_Area (thread_wa, 0);
        for  (i = 1;  i < Num_PThreads;  i ++)
          Initialize_Work_Area (thread_wa + i, i);
       }

   myRead = new_fragRecord ();

#if  SHOW_PROGRESS
Start_Time = clock ();
#endif

   if  (Contig_Mode)
       Next_Fragment_Index = 1;
     else
       Next_Fragment_Index = getLastElemFragStore (OldFragStore) + 1;



#if  USE_SOURCE_FIELD
Source_Log_File = File_Open ("ovl-srcinfo.log", "w");
#endif

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

   while (ReadFrags (Max_Hash_Strings))
     {
      int startIndex;

      if  (noOverlaps == 0)
          {
           GateKeeperStore  *curr_frag_store;
           GateKeeperStore  *hash_frag_store;
           int  highest_old_frag, lowest_old_frag;
           int  status;

                if  (Contig_Mode)
                    {
                     hash_frag_store
                         = loadFragStorePartial (BACtig_Store_Path,
                                                 First_Hash_Frag,
                                                 Last_Hash_Frag,
                                                 FRAG_S_INF | FRAG_S_SEQ | FRAG_S_QLT);
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
                                                 Last_Hash_Frag,
                                                 FRAG_S_INF | FRAG_S_SEQ | FRAG_S_QLT);
                     fprintf (stderr,
                              "loadFragStorePartial  first = " F_S64 "  last = " F_S64 "\n",
                              First_Hash_Frag, Last_Hash_Frag);
                     assert (0 < First_Hash_Frag
                               && First_Hash_Frag <= Last_Hash_Frag
                               && Last_Hash_Frag
                                    <= getLastElemFragStore (OldFragStore));
                    }
                HashFragStream = openFragStream (hash_frag_store, FRAG_S_INF | FRAG_S_SEQ | FRAG_S_QLT);
                resetFragStream (HashFragStream, First_Hash_Frag, Last_Hash_Frag);
                startIndex = First_Hash_Frag;


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



           lowest_old_frag = getFirstElemFragStore (OldFragStore);
           highest_old_frag = getLastElemFragStore (OldFragStore);
           if  (lowest_old_frag < Lo_Old_Frag)
               lowest_old_frag = Lo_Old_Frag;
           if  (highest_old_frag > Hi_Old_Frag)
               highest_old_frag = Hi_Old_Frag;
           if  (! Contig_Mode
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
                   Frag_Segment_Hi = Frag_Segment_Lo + Max_Frags_In_Memory_Store - 1;
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
                        (Frag_Store_Path, Frag_Segment_Lo, Frag_Segment_Hi, FRAG_S_INF | FRAG_S_SEQ | FRAG_S_QLT);
              fprintf (stderr,
                       "loadFragStorePartial  first = %d  last = %d\n",
                       Frag_Segment_Lo, Frag_Segment_Hi);
              assert (0 < Frag_Segment_Lo
                        && Frag_Segment_Lo <= Frag_Segment_Hi
                        && Frag_Segment_Hi
                             <= getLastElemFragStore (OldFragStore));

              for  (i = 0;  i < Num_PThreads;  i ++)
                {
                 old_stream_segment [i] = openFragStream (curr_frag_store, FRAG_S_INF | FRAG_S_SEQ | FRAG_S_QLT);
                 resetFragStream (old_stream_segment [i], Frag_Segment_Lo, Frag_Segment_Hi);
                }


              Now = time (NULL);
              fprintf (stderr, "### starting old fragments   %s\n", ctime (& Now));
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

              thread_wa [0] . stream_segment = old_stream_segment [0];
              Choose_And_Process_Stream_Segment (thread_wa);

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

              Now = time (NULL);
              fprintf (stderr, "### done old fragments   %s", ctime (& Now));
              for  (i = 0;  i < Num_PThreads;  i ++)
                closeFragStream (old_stream_segment [i]);
              closeGateKeeperStore (curr_frag_store);

              if  (IID_List == NULL)
                  lowest_old_frag += Max_Frags_In_Memory_Store;
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
                 
           closeFragStream (HashFragStream);
           closeGateKeeperStore (hash_frag_store);
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


//   Cleanup_Work_Area (driver_wa);
   Cleanup_Work_Area (thread_wa);
   for  (i = 1;  i < Num_PThreads;  i ++)
     Cleanup_Work_Area (thread_wa + i);
//   safe_free (driver_wa);
   safe_free (thread_wa);
   safe_free (thread_id);
   safe_free (new_stream_segment);
   safe_free (old_stream_segment);

#if  ANALYZE_HITS
fclose (High_Hits_File);
#endif

#if  USE_SOURCE_FIELD
fclose (Source_Log_File);
#endif

   closeFragStream (HashFragStream);

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
             && List [i] - List [lo] < Max_Frags_In_Memory_Store;
           i ++)
     ;

   return  i - 1;
  }



/******************************************************************************/

void  Cleanup_Work_Area
    (Work_Area_t * wa)

//  Free memory allocated in  (* wa) .

  {
   safe_free (wa -> String_Olap_Space);
   safe_free (wa -> Match_Node_Space);
   del_fragRecord (wa -> myRead);

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

      if  (Num_PThreads > 1)
          pthread_mutex_lock (& Fragment_Range_Mutex);

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

      if  (Num_PThreads > 1)
          pthread_mutex_unlock (& Fragment_Range_Mutex);

      if  (IID_List != NULL && hi > Frag_Segment_Hi)
          break;
      if  (hi > Frag_Segment_Hi)
          hi = Frag_Segment_Hi;
      if  (lo > hi)
          break;
      
      if  (Num_PThreads > 1)
          pthread_mutex_lock (& FragStore_Mutex);

      resetFragStream (WA -> stream_segment, lo, hi);

      if  (Num_PThreads > 1)
          pthread_mutex_unlock (& FragStore_Mutex);

      Process_Overlaps (WA -> stream_segment, WA);
     }

   if  (Num_PThreads > 1 && WA -> thread_id > 0)
       pthread_exit (ptr);

   return  ptr;
  }



/******************************************************************************/

/* function ReadFrags:
   Read up to maxFrags, but maybe a little less.

 */
static int  ReadFrags
    (int maxFrags)

  {
    Now = time (NULL);
    fprintf (stderr, "### Start batch #%d   %s", ++ Batch_Num, ctime (& Now));

    // In LSF mode there is no input file.  Instead, the hash table is built
    // from fragments from the fragment store.  In this case, we just
    // set values of global variables  First_Hash_Frag  and  Last_Hash_Frag
    // to indicate the range of fragments to get from that store.

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
