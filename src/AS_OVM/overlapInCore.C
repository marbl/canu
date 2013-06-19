
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

const char *mainid = "$Id: overlapInCore.C,v 1.12 2013-01-11 11:11:04 brianwalenz Exp $";

#include "overlapInCore.H"
#include "AS_UTL_decodeRange.H"


uint32 STRING_NUM_BITS       = 31;  //  MUST BE EXACTLY THIS
uint32 OFFSET_BITS           = 31;

uint64 TRUELY_ZERO           = 0;
uint64 TRUELY_ONE            = 1;

uint64 STRING_NUM_MASK       = (TRUELY_ONE << STRING_NUM_BITS) - 1;
uint64 OFFSET_MASK           = (TRUELY_ONE << OFFSET_BITS) - 1;

uint64 MAX_STRING_NUM        = STRING_NUM_MASK;



int64  Bad_Short_Window_Ct = 0;
//  The number of overlaps rejected because of too many errors in a small window

int64  Bad_Long_Window_Ct = 0;
//  The number of overlaps rejected because of too many errors in a long window

double  Branch_Match_Value = 0.0;
double  Branch_Error_Value = 0.0;
//  Scores of matches and mismatches in alignments.  Alignment
//  ends at maximum score.
//  THESE VALUES CAN BE SET ONLY AT RUN TIME (search for them below)

char  * Data = NULL;
//  Stores sequence data of fragments in hash table

char  * Quality_Data = NULL;
//  Stores quality data of fragments in hash table

size_t  Data_Len = 0;
bool  Doing_Partial_Overlaps = FALSE;
//  If set true by the G option (G for Granger)
//  then allow overlaps that do not extend to the end
//  of either read.

size_t  Extra_Data_Len;
//  Total length available for hash table string data,
//  including both regular strings and extra strings
//  added from kmer screening

uint64  Extra_Ref_Ct = 0;
String_Ref_t  * Extra_Ref_Space = NULL;
uint64  Extra_String_Ct = 0;
//  Number of extra strings of screen kmers added to hash table

uint64  Extra_String_Subcount = 0;
//  Number of kmers already added to last extra string in hash table

uint64  Frag_Olap_Limit = FRAG_OLAP_LIMIT;
//  Maximum number of overlaps for end of an old fragment against
//  a single hash table of frags, in each orientation

Check_Vector_t  * Hash_Check_Array = NULL;
//  Bit vector to eliminate impossible hash matches

uint64  Hash_String_Num_Offset = 1;
Hash_Bucket_t  * Hash_Table;
bool  Ignore_Clear_Range = FALSE;
//  If true will use entire read sequence, ignoring the
//  clear range values

uint64  Kmer_Hits_With_Olap_Ct = 0;
uint64  Kmer_Hits_Without_Olap_Ct = 0;
int32  Min_Olap_Len = 0;
uint64  Multi_Overlap_Ct = 0;
String_Ref_t  * Next_Ref = NULL;
uint64  String_Ct;
//  Number of fragments in the hash table

Hash_Frag_Info_t  * String_Info = NULL;
int64  * String_Start = NULL;
uint32  String_Start_Size = 0;
//  Number of available positions in  String_Start

bool  Unique_Olap_Per_Pair = TRUE;
//  If true will allow at most
//  one overlap output message per oriented fragment pair
//  Set true by  -u  command-line option; set false by  -m

size_t  Used_Data_Len = 0;
//  Number of bytes of Data currently occupied, including
//  regular strings and extra kmer screen strings

bool  Use_Hopeless_Check = TRUE;
//  Determines whether check for absence of kmer matches
//  at the end of a read is used to abort the overlap before
//  the extension from a single kmer match is attempted.

bool  Use_Window_Filter = FALSE;
//  Determines whether check for a window containing too many
//  errors is used to disqualify overlaps.

int32  Read_Edit_Match_Limit [AS_READ_MAX_NORMAL_LEN] = {0};
//  This array [e] is the minimum value of  Edit_Array [e] [d]
//  to be worth pursuing in edit-distance computations between reads
//  (only MAX_ERRORS needed)

int32  Guide_Edit_Match_Limit [AS_READ_MAX_NORMAL_LEN] = {0};
//  This array [e] is the minimum value of  Edit_Array [e] [d]
//  to be worth pursuing in edit-distance computations between guides
//  (only MAX_ERRORS needed)

int32  Read_Error_Bound [AS_READ_MAX_NORMAL_LEN + 1];
//  This array [i]  is the maximum number of errors allowed
//  in a match between reads of length  i , which is
//  i * AS_OVL_ERROR_RATE .

int32  Guide_Error_Bound [AS_READ_MAX_NORMAL_LEN + 1];
//  This array [i]  is the maximum number of errors allowed
//  in a match between guides of length  i , which is
//  i * AS_OVL_ERROR_RATE .

double  Branch_Cost [AS_READ_MAX_NORMAL_LEN + 1];
//  Branch_Cost [i]  is the "goodness" of matching i characters
//  after a single error in determining branch points.

int32  Bit_Equivalent [256] = {0};
//  Table to convert characters to 2-bit integer code

int32  Char_Is_Bad [256] = {0};
//  Table to check if character is not a, c, g or t.

uint64  Hash_Entries = 0;

uint64  Total_Overlaps = 0;
uint64  Contained_Overlap_Ct = 0;
uint64  Dovetail_Overlap_Ct = 0;

uint32 minLibToHash = 0;
uint32 maxLibToHash = 0;
uint32 minLibToRef  = 0;
uint32 maxLibToRef  = 0;

uint64  Kmer_Len = 0;
uint64  HSF1     = 666;
uint64  HSF2     = 666;
uint64  SV1      = 666;
uint64  SV2      = 666;
uint64  SV3      = 666;

uint32  Hash_Mask_Bits            = 22;
double  Max_Hash_Load             = 0.6;
uint32  Max_Hash_Strings          = 100000000 / 800;
uint64  Max_Hash_Data_Len         = 100000000;

uint32  Max_Reads_Per_Batch       = 0;  //  The number of reads loaded in a single batch.
uint32  Max_Reads_Per_Thread      = 0;  //  The number of reads processed per thread.

AS_IID  Last_Hash_Frag_Read;
AS_IID  Lo_Hash_Frag = 0;
AS_IID  Hi_Hash_Frag = AS_IID_MAX;
AS_IID  Lo_Old_Frag  = 0;
AS_IID  Hi_Old_Frag  = AS_IID_MAX;
uint32  Num_PThreads = 4;

char  Sequence_Buffer [2 * AS_READ_MAX_NORMAL_LEN];
char  Quality_Buffer [2 * AS_READ_MAX_NORMAL_LEN];

FILE  * Kmer_Skip_File = NULL;

BinaryOverlapFile  *Out_BOF = NULL;

gkStore  *OldFragStore;
char  * Frag_Store_Path;

pthread_mutex_t  FragStore_Mutex;
pthread_mutex_t  Write_Proto_Mutex;

AS_IID      First_Hash_Frag = 0;
AS_IID      Last_Hash_Frag  = 0;
gkFragment  myRead;

AS_IID  Frag_Segment_Lo;
AS_IID  Frag_Segment_Hi;



//  Find all overlaps between frags  Frag_Segment_Lo .. Frag_Segment_Hi
//  in the stream in  ptr  and the frags in the hash table.
//
static
void *Choose_And_Process_Stream_Segment(void *ptr) {
  Work_Area_t  *WA = (Work_Area_t *) (ptr);
  int           allDone = 0;

  fprintf(stderr, "Choose_And_Process_Stream_Segment()-- tid %d\n", WA->thread_id);

  while  (allDone == 0) {
    pthread_mutex_lock (& FragStore_Mutex);

    AS_IID lo = Frag_Segment_Lo;
    Frag_Segment_Lo += Max_Reads_Per_Thread;
    AS_IID hi = Frag_Segment_Lo - 1;

    //  This block doesn't need to be in a mutex, but it's so quick
    //  we just do it rather than exiting and entering a mutex
    //  again.
    //
    if  (hi > Frag_Segment_Hi)
      hi = Frag_Segment_Hi;
    if  (lo > hi)
      allDone = 1;

    //  This definitely needs to be mutex'd.
    if (allDone == 0)
      WA->stream_segment->reset (lo, hi);

    pthread_mutex_unlock (& FragStore_Mutex);

    if (allDone == 0)
      Process_Overlaps (WA -> stream_segment, WA);
  }

  fprintf(stderr, "Choose_And_Process_Stream_Segment()-- tid %d returns\n", WA->thread_id);

  return(ptr);
}



//  Allocate memory for  (* WA)  and set initial values.
//  Set  thread_id  field to  id .
void
Initialize_Work_Area(Work_Area_t * WA, int id) {
  uint64  allocated = 0;

  WA -> Left_Delta  = (int *)safe_malloc (MAX_ERRORS * sizeof (int));
  WA -> Right_Delta = (int *)safe_malloc (MAX_ERRORS * sizeof (int));

  WA -> Delta_Stack = (int *)safe_malloc (MAX_ERRORS * sizeof (int));

  allocated += 3 * MAX_ERRORS * sizeof(int);

  WA -> Edit_Space = (int *)safe_malloc ((MAX_ERRORS + 4) * MAX_ERRORS * sizeof (int));
  WA -> Edit_Array = (int **)safe_malloc (MAX_ERRORS * sizeof (int *));

  allocated += (MAX_ERRORS + 4) * MAX_ERRORS * sizeof (int) + MAX_ERRORS * sizeof (int *);

  WA -> String_Olap_Size = INIT_STRING_OLAP_SIZE;
  WA -> String_Olap_Space = (String_Olap_t *) safe_malloc (WA -> String_Olap_Size * sizeof (String_Olap_t));
  WA -> Match_Node_Size = INIT_MATCH_NODE_SIZE;
  WA -> Match_Node_Space = (Match_Node_t *) safe_malloc (WA -> Match_Node_Size * sizeof (Match_Node_t));

  allocated += WA -> String_Olap_Size * sizeof (String_Olap_t);
  allocated += WA -> Match_Node_Size * sizeof (Match_Node_t);

  int32 Offset = 2;
  int32 Del = 6;
  for  (int32 i=0;  i<MAX_ERRORS;  i++) {
    WA -> Edit_Array [i] = WA -> Edit_Space + Offset;
    Offset += Del;
    Del += 2;
  }

  WA -> status = 0;
  WA -> thread_id = id;

  //  OVSoverlap is 16 bytes, so 1MB of data would store 65536
  //  overlaps.
  //
  WA->overlapsLen = 0;
  WA->overlapsMax = 1024 * 1024 / sizeof(OVSoverlap);
  WA->overlaps    = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * WA->overlapsMax);

  allocated += sizeof(OVSoverlap) * WA->overlapsMax;

  fprintf(stderr, "Initialize_Work_Area:  MAX_ERRORS=%d  allocated "F_U64"MB\n", MAX_ERRORS, allocated >> 20);
}





static
bool
ReadFrags(AS_IID maxFrags) {

  if  (First_Hash_Frag == 0)
    First_Hash_Frag = Lo_Hash_Frag;
  else
    First_Hash_Frag = Last_Hash_Frag + 1;

  if (First_Hash_Frag > Hi_Hash_Frag)
    return(false);

  Last_Hash_Frag = First_Hash_Frag + maxFrags - 1;

  if (Hi_Hash_Frag < Last_Hash_Frag)
    Last_Hash_Frag = Hi_Hash_Frag;

  return(true);
}




int
OverlapDriver(void) {
  pthread_attr_t  attr;
  pthread_t  * thread_id;
  gkStream **new_stream_segment;
  gkStream **old_stream_segment;
  Work_Area_t  * thread_wa;

  thread_id = (pthread_t *) safe_calloc (Num_PThreads, sizeof (pthread_t));

  new_stream_segment = (gkStream **) safe_calloc (Num_PThreads, sizeof (gkStream *));
  old_stream_segment = (gkStream **) safe_calloc (Num_PThreads, sizeof (gkStream *));
  thread_wa = (Work_Area_t *) safe_calloc (Num_PThreads, sizeof (Work_Area_t));

  for  (uint32 i = 0;  i < Num_PThreads;  i ++)
    old_stream_segment [i] = new gkStream (OldFragStore, 0, 0, GKFRAGMENT_QLT);

  pthread_attr_init (& attr);
  pthread_attr_setstacksize (& attr, THREAD_STACKSIZE);
  pthread_mutex_init (& FragStore_Mutex, NULL);
  pthread_mutex_init (& Write_Proto_Mutex, NULL);

  Initialize_Work_Area (thread_wa, 0);
  for  (uint32 i = 1;  i < Num_PThreads;  i ++)
    Initialize_Work_Area (thread_wa + i, i);

  {
    AS_IID  id;

    if  (Lo_Hash_Frag < 1)
      Lo_Hash_Frag = 1;

    id = OldFragStore->gkStore_getNumFragments ();

    if  (id < Hi_Hash_Frag)
      Hi_Hash_Frag = id;
  }

  while (ReadFrags (Max_Hash_Strings)) {
    gkStore  *curr_frag_store;
    gkStore  *hash_frag_store;

    hash_frag_store = new gkStore(Frag_Store_Path, FALSE, FALSE);
    hash_frag_store->gkStore_load(First_Hash_Frag, Last_Hash_Frag, GKFRAGMENT_QLT);
    assert (0 < First_Hash_Frag
            && First_Hash_Frag <= Last_Hash_Frag
            && Last_Hash_Frag  <= OldFragStore->gkStore_getNumFragments ());

    fprintf(stderr, "Build_Hash_Index from "F_IID" to "F_IID"\n", First_Hash_Frag, Last_Hash_Frag);

    gkStream *hashStream = new gkStream (hash_frag_store, First_Hash_Frag, Last_Hash_Frag, GKFRAGMENT_QLT);
    Build_Hash_Index (hashStream, First_Hash_Frag, &myRead);
    delete hashStream;

    if (Last_Hash_Frag_Read < Last_Hash_Frag)
      //  Didn't read all frags.
      Last_Hash_Frag = Last_Hash_Frag_Read;

    fprintf(stderr, "Index built.\n");

    AS_IID lowest_old_frag  = 1;
    AS_IID highest_old_frag = OldFragStore->gkStore_getNumFragments ();

    if  (lowest_old_frag < Lo_Old_Frag)
      lowest_old_frag = Lo_Old_Frag;
    if  (highest_old_frag > Hi_Old_Frag)
      highest_old_frag = Hi_Old_Frag;
    if  (highest_old_frag > Last_Hash_Frag)
      highest_old_frag = Last_Hash_Frag;

    while  (lowest_old_frag <= highest_old_frag) {
      Frag_Segment_Lo = lowest_old_frag;
      Frag_Segment_Hi = Frag_Segment_Lo + Max_Reads_Per_Batch - 1;
      if  (Frag_Segment_Hi > highest_old_frag)
        Frag_Segment_Hi = highest_old_frag;

      fprintf(stderr, "Starting "F_U32" "F_U32"\n", Frag_Segment_Lo, Frag_Segment_Hi);

      curr_frag_store = new gkStore(Frag_Store_Path, FALSE, FALSE);
      curr_frag_store->gkStore_load(Frag_Segment_Lo, Frag_Segment_Hi, GKFRAGMENT_QLT);
      assert(0 < Frag_Segment_Lo);
      assert(Frag_Segment_Lo <= Frag_Segment_Hi);
      assert(Frag_Segment_Hi <= OldFragStore->gkStore_getNumFragments ());

      for  (uint32 i = 0;  i < Num_PThreads;  i ++)
        old_stream_segment [i] = new gkStream (curr_frag_store, Frag_Segment_Lo, Frag_Segment_Hi, GKFRAGMENT_QLT);

      for  (uint32 i = 1;  i < Num_PThreads;  i ++) {
        thread_wa [i] . stream_segment = old_stream_segment [i];
        int status = pthread_create (thread_id + i, & attr, Choose_And_Process_Stream_Segment, thread_wa + i);
        if  (status != 0)
          fprintf (stderr, "pthread_create error:  %s\n", strerror(status)), exit(1);
      }

      thread_wa [0] . stream_segment = old_stream_segment [0];
      Choose_And_Process_Stream_Segment (thread_wa);

      for  (uint32 i = 1;  i < Num_PThreads;  i ++) {
        int status = pthread_join  (thread_id [i], NULL);
        if (status != 0)
          fprintf(stderr, "pthread_join error: %s\n", strerror(status)), exit(1);
      }

      for  (uint32 i = 0;  i < Num_PThreads;  i ++)
        delete old_stream_segment [i];

      delete curr_frag_store;

      lowest_old_frag += Max_Reads_Per_Batch;
    }

    delete hash_frag_store;
  }

  for  (uint32 i=0;  i<Num_PThreads;  i++) {
    safe_free (thread_wa[i].String_Olap_Space);
    safe_free (thread_wa[i].Match_Node_Space);
    safe_free (thread_wa[i].overlaps);
  }

  safe_free (thread_wa);
  safe_free (thread_id);

  safe_free (new_stream_segment);
  safe_free (old_stream_segment);

  return  0;
}










//  Return the smallest  n >= Start  s.t.
//    prob [>= e  errors in  n  binomial trials (p = error prob)]
//          > Limit

static
int
Binomial_Bound (int e, double p, int Start, double Limit) {
  double  Normal_Z, Mu_Power, Factorial, Poisson_Coeff;
  double  q, Sum, P_Power, Q_Power, X;
  int  k, n, Bin_Coeff, Ct;

  q = 1.0 - p;
  if  (Start < e)
    Start = e;

  for  (n = Start;  n < AS_READ_MAX_NORMAL_LEN;  n ++) {
    if  (n <= 35) {
      Sum = 0.0;
      Bin_Coeff = 1;
      Ct = 0;
      P_Power = 1.0;
      Q_Power = pow (q, n);

      for  (k = 0;  k < e && 1.0 - Sum > Limit;  k ++) {
        X = Bin_Coeff * P_Power * Q_Power;
        Sum += X;
        Bin_Coeff *= n - Ct;
        Bin_Coeff /= ++ Ct;
        P_Power *= p;
        Q_Power /= q;
      }
      if  (1.0 - Sum > Limit)
        return  n;
    } else {
      Normal_Z = (e - 0.5 - n * p) / sqrt (n * p * q);
      if  (Normal_Z <= NORMAL_DISTRIB_THOLD)
        return  n;
#undef COMPUTE_IN_LOG_SPACE
#ifndef COMPUTE_IN_LOG_SPACE
      Sum = 0.0;
      Mu_Power = 1.0;
      Factorial = 1.0;
      Poisson_Coeff = exp (- n * p);
      for  (k = 0;  k < e;  k ++) {
        Sum += Mu_Power * Poisson_Coeff / Factorial;
        Mu_Power *= n * p;
        Factorial *= k + 1;
      }
#else
      Sum = 0.0;
      Mu_Power = 0.0;
      Factorial = 0.0;
      Poisson_Coeff = - n * p;
      for  (k = 0;  k < e;  k ++) {
	      Sum += exp(Mu_Power + Poisson_Coeff - Factorial);
        Mu_Power += log(n * p);
        Factorial = lgamma(k + 1);
      }
#endif
      if  (1.0 - Sum > Limit)
        return  n;
    }
  }

  return  AS_READ_MAX_NORMAL_LEN;
}


static
void
Initialize_Globals (void) {

  for  (int i = 0;  i <= ERRORS_FOR_FREE;  i ++)
    Read_Edit_Match_Limit [i] = 0;

  int Start = 1;

  for  (int e = ERRORS_FOR_FREE + 1;  e < MAX_ERRORS;  e ++) {
    Start = Binomial_Bound (e - ERRORS_FOR_FREE, AS_OVL_ERROR_RATE, Start, EDIT_DIST_PROB_BOUND);
    Read_Edit_Match_Limit [e] = Start - 1;
    assert (Read_Edit_Match_Limit [e] >= Read_Edit_Match_Limit [e - 1]);
  }

  for  (int i = 0;  i <= ERRORS_FOR_FREE;  i ++)
    Guide_Edit_Match_Limit [i] = 0;

  Start = 1;

  for  (int e = ERRORS_FOR_FREE + 1;  e < MAX_ERRORS;  e ++) {
    Start = Binomial_Bound (e - ERRORS_FOR_FREE, AS_OVL_ERROR_RATE, Start, EDIT_DIST_PROB_BOUND);
    Guide_Edit_Match_Limit [e] = Start - 1;
    assert (Guide_Edit_Match_Limit [e] >= Guide_Edit_Match_Limit [e - 1]);
  }

  for  (int i = 0;  i <= AS_READ_MAX_NORMAL_LEN;  i ++)
    Read_Error_Bound [i] = (int) (i * AS_OVL_ERROR_RATE + 0.0000000000001);

  for  (int i = 0;  i <= AS_READ_MAX_NORMAL_LEN;  i ++)
    Guide_Error_Bound [i] = (int) (i * AS_OVL_ERROR_RATE + 0.0000000000001);

  for  (int i = 0;  i <= AS_READ_MAX_NORMAL_LEN;  i ++)
    Branch_Cost [i] = i * Branch_Match_Value + Branch_Error_Value;

  Bit_Equivalent ['a'] = Bit_Equivalent ['A'] = 0;
  Bit_Equivalent ['c'] = Bit_Equivalent ['C'] = 1;
  Bit_Equivalent ['g'] = Bit_Equivalent ['G'] = 2;
  Bit_Equivalent ['t'] = Bit_Equivalent ['T'] = 3;

  for  (int i = 0;  i < 256;  i ++) {
    char  ch = tolower ((char) i);

    if  (ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
      Char_Is_Bad [i] = 0;
    else
      Char_Is_Bad [i] = 1;
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "HASH_TABLE_SIZE         "F_U32"\n",     HASH_TABLE_SIZE);
  fprintf(stderr, "sizeof(Hash_Bucket_t)   "F_SIZE_T"\n",  sizeof(Hash_Bucket_t));
  fprintf(stderr, "hash table size:        "F_SIZE_T" MB\n",  (HASH_TABLE_SIZE * sizeof(Hash_Bucket_t)) >> 20);
  fprintf(stderr, "\n");

  Hash_Table = (Hash_Bucket_t *) safe_malloc (HASH_TABLE_SIZE * sizeof (Hash_Bucket_t));

  fprintf(stderr, "check  "F_SIZE_T" MB\n", (HASH_TABLE_SIZE * sizeof (Check_Vector_t) >> 20));
  fprintf(stderr, "info   "F_SIZE_T" MB\n", (Max_Hash_Strings * sizeof (Hash_Frag_Info_t) >> 20));
  fprintf(stderr, "start  "F_SIZE_T" MB\n", (Max_Hash_Strings * sizeof (int64) >> 20));
  fprintf(stderr, "\n");

  Hash_Check_Array = (Check_Vector_t *) safe_malloc (HASH_TABLE_SIZE * sizeof (Check_Vector_t));
  String_Info = (Hash_Frag_Info_t *) safe_calloc (Max_Hash_Strings, sizeof (Hash_Frag_Info_t));
  String_Start = (int64 *) safe_calloc (Max_Hash_Strings, sizeof (int64));

  String_Start_Size = Max_Hash_Strings;
}





int
main(int argc, char **argv) {
  char  bolfile_name[FILENAME_MAX] = {0};
  char  Outfile_Name[FILENAME_MAX] = {0};
  int  illegal;
  char  * p;

  argc = AS_configure(argc, argv);
  Min_Olap_Len = AS_OVERLAP_MIN_LEN; // set after configure

  int err=0;
  int arg=1;
  while (arg < argc) {
    if (strcmp(argv[arg], "-G") == 0) {
      Doing_Partial_Overlaps = TRUE;
    } else if (strcmp(argv[arg], "-h") == 0) {
      AS_UTL_decodeRange(argv[++arg], Lo_Hash_Frag, Hi_Hash_Frag);

    } else if (strcmp(argv[arg], "-H") == 0) {
      AS_UTL_decodeRange(argv[++arg], minLibToHash, maxLibToHash);

    } else if (strcmp(argv[arg], "-R") == 0) {
      AS_UTL_decodeRange(argv[++arg], minLibToRef, maxLibToRef);

    } else if (strcmp(argv[arg], "-k") == 0) {
      arg++;
      if ((isdigit(argv[arg][0]) && (argv[arg][1] == 0)) ||
          (isdigit(argv[arg][0]) && isdigit(argv[arg][1]) && (argv[arg][2] == 0))) {
        Kmer_Len = strtoull(argv[arg], NULL, 10);
      } else {
        errno = 0;
        Kmer_Skip_File = fopen(argv[arg], "r");
        if (errno)
          fprintf(stderr, "ERROR: Failed to open -k '%s': %s\n", argv[arg], strerror(errno)), exit(1);
      }

    } else if (strcmp(argv[arg], "-l") == 0) {
      Frag_Olap_Limit = strtol(argv[++arg], NULL, 10);
      if  (Frag_Olap_Limit < 1)
        Frag_Olap_Limit = INT_MAX;

    } else if (strcmp(argv[arg], "-m") == 0) {
      Unique_Olap_Per_Pair = FALSE;

    } else if (strcmp(argv[arg], "--hashbits") == 0) {
      Hash_Mask_Bits = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--hashstrings") == 0) {
      Max_Hash_Strings = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--hashdatalen") == 0) {
      Max_Hash_Data_Len = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--hashload") == 0) {
      Max_Hash_Load = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "--maxreadlen") == 0) {
      //  Quite the gross way to do this, but simple.
      uint32 desired = strtoul(argv[++arg], NULL, 10);
      OFFSET_BITS = 1;
      while (((uint32)1 << OFFSET_BITS) < desired)
        OFFSET_BITS++;

      STRING_NUM_BITS       = 30 - OFFSET_BITS;

      STRING_NUM_MASK       = (1 << STRING_NUM_BITS) - 1;
      OFFSET_MASK           = (1 << OFFSET_BITS) - 1;

      MAX_STRING_NUM        = STRING_NUM_MASK;

    } else if (strcmp(argv[arg], "--readsperbatch") == 0) {
      Max_Reads_Per_Batch = strtoul(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--readsperthread") == 0) {
      Max_Reads_Per_Thread = strtoul(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-o") == 0) {
      strcpy(Outfile_Name, argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      AS_UTL_decodeRange(argv[++arg], Lo_Old_Frag, Hi_Old_Frag);

    } else if (strcmp(argv[arg], "-t") == 0) {
      Num_PThreads = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-u") == 0) {
      Unique_Olap_Per_Pair = TRUE;

    } else if (strcmp(argv[arg], "-v") == 0) {
      Min_Olap_Len = (int) strtol (argv[++arg], & p, 10);

    } else if (strcmp(argv[arg], "-w") == 0) {
      Use_Window_Filter = TRUE;

    } else if (strcmp(argv[arg], "-x") == 0) {
      Ignore_Clear_Range = TRUE;

    } else if (strcmp(argv[arg], "-z") == 0) {
      Use_Hopeless_Check = FALSE;

    } else {
      if (Frag_Store_Path == NULL) {
        Frag_Store_Path = argv[arg];
      } else {
        fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
        err++;
      }
    }
    arg++;
  }

  //  Fix up some flags if we're allowing high error rates.
  //
  if (AS_OVL_ERROR_RATE > 0.06) {
    if (Use_Window_Filter)
      fprintf(stderr, "High error rates requested -- window-filter turned off despite -w flag!\n");
    Use_Window_Filter  = FALSE;
    Use_Hopeless_Check = FALSE;
  }

  if (Max_Hash_Strings == 0)
    fprintf(stderr, "* No memory model supplied; -M needed!\n"), err++;

  if (Kmer_Len == 0)
    fprintf(stderr, "* No kmer length supplied; -k needed!\n"), err++;

  if (Max_Hash_Strings > MAX_STRING_NUM)
    fprintf(stderr, "Too many strings (--hashstrings), must be less than "F_U64"\n", MAX_STRING_NUM), err++;

  if (Outfile_Name[0] == 0)
    fprintf (stderr, "ERROR:  No output file name specified\n"), err++;

  if ((err) || (Frag_Store_Path == NULL)) {
    fprintf(stderr, "USAGE:  %s [options] <gkpStorePath>\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "-b <fn>     in contig mode, specify the output file\n");
    fprintf(stderr, "-c          contig mode.  Use 2 frag stores.  First is\n");
    fprintf(stderr, "            for reads; second is for contigs\n");
    fprintf(stderr, "-G          do partial overlaps\n");
    fprintf(stderr, "-h <range>  to specify fragments to put in hash table\n");
    fprintf(stderr, "            Implies LSF mode (no changes to frag store)\n");
    fprintf(stderr, "-I          designate a file of frag iids to limit olaps to\n");
    fprintf(stderr, "            (Contig mode only)\n");
    fprintf(stderr, "-k          if one or two digits, the length of a kmer, otherwise\n");
    fprintf(stderr, "            the filename containing a list of kmers to ignore in\n");
    fprintf(stderr, "            the hash table\n");
    fprintf(stderr, "-l          specify the maximum number of overlaps per\n");
    fprintf(stderr, "            fragment-end per batch of fragments.\n");
    fprintf(stderr, "-m          allow multiple overlaps per oriented fragment pair\n");
    fprintf(stderr, "-M          specify memory size.  Valid values are '8GB', '4GB',\n");
    fprintf(stderr, "            '2GB', '1GB', '256MB'.  (Not for Contig mode)\n");
    fprintf(stderr, "-o          specify output file name\n");
    fprintf(stderr, "-P          write protoIO output (if not -G)\n");
    fprintf(stderr, "-r <range>  specify old fragments to overlap\n");
    fprintf(stderr, "-s          ignore screen information with fragments\n");
    fprintf(stderr, "-t <n>      use <n> parallel threads\n");
    fprintf(stderr, "-u          allow only 1 overlap per oriented fragment pair\n");
    fprintf(stderr, "-v <n>      only output overlaps of <n> or more bases\n");
    fprintf(stderr, "-w          filter out overlaps with too many errors in a window\n");
    fprintf(stderr, "-x          ignore the clear ranges on reads and use the \n");
    fprintf(stderr, "            full sequence\n");
    fprintf(stderr, "-z          skip the hopeless check\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--hashbits n       Use n bits for the hash mask.\n");
    fprintf(stderr, "--hashstrings n    Load at most n strings into the hash table at one time.\n");
    fprintf(stderr, "--hashdatalen n    Load at most n bytes into the hash table at one time.\n");
    fprintf(stderr, "--hashload f       Load to at most 0.0 < f < 1.0 capacity (default 0.7).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--maxreadlen n     For batches with all short reads, pack bits differently to\n");
    fprintf(stderr, "                   process more reads per batch.\n");
    fprintf(stderr, "                     all reads must be shorter than n\n");
    fprintf(stderr, "                     --hashstrings limited to 2^(30-m)\n");
    fprintf(stderr, "                   Common values:\n");
    fprintf(stderr, "                     maxreadlen 2048 -> hashstrings  524288 (default)\n");
    fprintf(stderr, "                     maxreadlen  512 -> hashstrings 2097152\n");
    fprintf(stderr, "                     maxreadlen  128 -> hashstrings 8388608\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--readsperbatch n  Force batch size to n.\n");
    fprintf(stderr, "--readsperthread n Force each thread to process n reads.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  assert(NULL == Out_BOF);

  Out_BOF    = AS_OVS_createBinaryOverlapFile(Outfile_Name, FALSE);

  //  Adjust the number of reads to load into memory at once (for processing, not the hash table),

  if (Max_Reads_Per_Batch == 0)
    Max_Reads_Per_Batch = 100000;

  if (Max_Hash_Strings < Max_Reads_Per_Batch)
    Max_Reads_Per_Batch = Max_Hash_Strings;

  //  Adjust the number of reads processed per thread.  Default to having four blocks per thread,
  //  but make sure that (a) all threads have work to do, and (b) batches are not minuscule.

  if (Max_Reads_Per_Thread == 0)
    Max_Reads_Per_Thread = Max_Reads_Per_Batch / (4 * Num_PThreads);

  if (Max_Reads_Per_Thread * Num_PThreads > Max_Reads_Per_Batch)
    Max_Reads_Per_Thread = Max_Reads_Per_Batch / Num_PThreads + 1;

  if (Max_Reads_Per_Thread < 10)
    Max_Reads_Per_Thread = 10;

  //  We know enough now to set the hash function variables, and some other random variables.

  HSF1 = Kmer_Len - (Hash_Mask_Bits / 2);
  HSF2 = 2 * Kmer_Len - Hash_Mask_Bits;
  SV1  = HSF1 + 2;
  SV2  = (HSF1 + HSF2) / 2;
  SV3  = HSF2 - 2;

  Branch_Match_Value = (Doing_Partial_Overlaps) ? PARTIAL_BRANCH_MATCH_VAL : DEFAULT_BRANCH_MATCH_VAL;
  Branch_Error_Value = Branch_Match_Value - 1.0;

  fprintf(stderr, "\n");
  fprintf(stderr, "STRING_NUM_BITS       "F_U32"\n", STRING_NUM_BITS);
  fprintf(stderr, "OFFSET_BITS           "F_U32"\n", OFFSET_BITS);
  fprintf(stderr, "STRING_NUM_MASK       "F_U64"\n", STRING_NUM_MASK);
  fprintf(stderr, "OFFSET_MASK           "F_U64"\n", OFFSET_MASK);
  fprintf(stderr, "MAX_STRING_NUM        "F_U64"\n", MAX_STRING_NUM);
  fprintf(stderr, "\n");
  fprintf(stderr, "Hash_Mask_Bits        "F_U32"\n", Hash_Mask_Bits);
  fprintf(stderr, "Max_Hash_Strings      "F_U32"\n", Max_Hash_Strings);
  fprintf(stderr, "Max_Hash_Data_Len     "F_U64"\n", Max_Hash_Data_Len);
  fprintf(stderr, "Max_Hash_Load         %f\n", Max_Hash_Load);
  fprintf(stderr, "Kmer Length           %d\n", (int)Kmer_Len);
  fprintf(stderr, "Min Overlap Length    %d\n", Min_Olap_Len);
  fprintf(stderr, "MAX_ERRORS            %d\n", MAX_ERRORS);
  fprintf(stderr, "ERRORS_FOR_FREE       %d\n", ERRORS_FOR_FREE);
  fprintf(stderr, "\n");
  fprintf(stderr, "Num_PThreads          "F_U32"\n", Num_PThreads);
  fprintf(stderr, "Max_Reads_Per_Batch   "F_U32"\n", Max_Reads_Per_Batch);
  fprintf(stderr, "Max_Reads_Per_Thread  "F_U32"\n", Max_Reads_Per_Thread);

  assert (8 * sizeof (uint64) > 2 * Kmer_Len);

  Initialize_Globals ();

  OldFragStore = new gkStore(Frag_Store_Path, FALSE, FALSE);

  /****************************************/
  OverlapDriver();
  /****************************************/

  fprintf (stderr, " Kmer hits without olaps = "F_S64"\n", Kmer_Hits_Without_Olap_Ct);
  fprintf (stderr, "    Kmer hits with olaps = "F_S64"\n", Kmer_Hits_With_Olap_Ct);
  fprintf (stderr, "  Multiple overlaps/pair = "F_S64"\n", Multi_Overlap_Ct);
  fprintf (stderr, " Total overlaps produced = "F_S64"\n", Total_Overlaps);
  fprintf (stderr, "      Contained overlaps = "F_S64"\n", Contained_Overlap_Ct);
  fprintf (stderr, "       Dovetail overlaps = "F_S64"\n", Dovetail_Overlap_Ct);
  fprintf (stderr, "Rejected by short window = "F_S64"\n", Bad_Short_Window_Ct);
  fprintf (stderr, " Rejected by long window = "F_S64"\n", Bad_Long_Window_Ct);

  delete OldFragStore;

  AS_OVS_closeBinaryOverlapFile(Out_BOF);

  return(0);
}
