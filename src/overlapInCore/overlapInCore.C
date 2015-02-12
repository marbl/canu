
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

const char *mainid = "$Id$";

#include "overlapInCore.H"
#include "AS_UTL_decodeRange.H"

oicParameters  G;


uint32 STRING_NUM_BITS       = 31;  //  MUST BE EXACTLY THIS
uint32 OFFSET_BITS           = 31;

uint64 TRUELY_ONE            = (uint64)1;
uint64 TRUELY_ZERO           = (uint64)0;

uint64 STRING_NUM_MASK       = (TRUELY_ONE << STRING_NUM_BITS) - 1;
uint64 OFFSET_MASK           = (TRUELY_ONE << OFFSET_BITS) - 1;

uint64 MAX_STRING_NUM        = STRING_NUM_MASK;



int64  Bad_Short_Window_Ct = 0;
//  The number of overlaps rejected because of too many errors in a small window

int64  Bad_Long_Window_Ct = 0;
//  The number of overlaps rejected because of too many errors in a long window


//  Stores sequence and quality data of fragments in hash table
char   *basesData = NULL;
char   *qualsData = NULL;
size_t  Data_Len = 0;

String_Ref_t  *nextRef = NULL;

size_t  Extra_Data_Len;
//  Total length available for hash table string data,
//  including both regular strings and extra strings
//  added from kmer screening

uint64         Max_Extra_Ref_Space = 0;  //  allocated amount
uint64         Extra_Ref_Ct = 0;         //  used amount
String_Ref_t  *Extra_Ref_Space = NULL;
uint64         Extra_String_Ct = 0;
//  Number of extra strings of screen kmers added to hash table

uint64  Extra_String_Subcount = 0;
//  Number of kmers already added to last extra string in hash table

Check_Vector_t  * Hash_Check_Array = NULL;
//  Bit vector to eliminate impossible hash matches

uint64  Hash_String_Num_Offset = 1;
Hash_Bucket_t  * Hash_Table;

uint64  Kmer_Hits_With_Olap_Ct = 0;
uint64  Kmer_Hits_Without_Olap_Ct = 0;
uint64  Multi_Overlap_Ct = 0;

uint64  String_Ct;
//  Number of fragments in the hash table

Hash_Frag_Info_t  * String_Info = NULL;
int64  * String_Start = NULL;
uint32  String_Start_Size = 0;
//  Number of available positions in  String_Start

size_t  Used_Data_Len = 0;
//  Number of bytes of Data currently occupied, including
//  regular strings and extra kmer screen strings

int32  Bit_Equivalent[256] = {0};
//  Table to convert characters to 2-bit integer code

int32  Char_Is_Bad[256] = {0};
//  Table to check if character is not a, c, g or t.

uint64  Hash_Entries = 0;

uint64  Total_Overlaps = 0;
uint64  Contained_Overlap_Ct = 0;
uint64  Dovetail_Overlap_Ct = 0;

uint64  HSF1     = 666;
uint64  HSF2     = 666;
uint64  SV1      = 666;
uint64  SV2      = 666;
uint64  SV3      = 666;

pthread_mutex_t  Write_Proto_Mutex;

ovFile  *Out_BOF = NULL;



//  Allocate memory for  (* WA)  and set initial values.
//  Set  thread_id  field to  id .
void
Initialize_Work_Area(Work_Area_t *WA, int id) {
  uint64  allocated = 0;

  WA->String_Olap_Size  = INIT_STRING_OLAP_SIZE;
  WA->String_Olap_Space = new String_Olap_t [WA->String_Olap_Size];

  WA->Match_Node_Size  = INIT_MATCH_NODE_SIZE;
  WA->Match_Node_Space = new Match_Node_t [WA->Match_Node_Size];

  allocated += WA->String_Olap_Size * sizeof (String_Olap_t);
  allocated += WA->Match_Node_Size  * sizeof (Match_Node_t);

  WA->status     = 0;
  WA->thread_id  = id;

  WA->overlapsLen = 0;
  WA->overlapsMax = 1024 * 1024 / sizeof(ovsOverlap);
  WA->overlaps    = new ovsOverlap [WA->overlapsMax];

  allocated += sizeof(ovsOverlap) * WA->overlapsMax;

  WA->editDist = new prefixEditDistance(G.Doing_Partial_Overlaps);
}


void
Delete_Work_Area(Work_Area_t *WA) {
  delete    WA->editDist;
  delete [] WA->String_Olap_Space;
  delete [] WA->Match_Node_Space;
  delete [] WA->overlaps;
}




int
OverlapDriver(void) {

  pthread_t      *thread_id = new pthread_t   [G.Num_PThreads];
  Work_Area_t    *thread_wa = new Work_Area_t [G.Num_PThreads];

  gkStore        *gkpStore  = new gkStore(G.Frag_Store_Path);

  pthread_attr_t  attr;

  pthread_attr_init(&attr);
  pthread_attr_setstacksize(&attr, THREAD_STACKSIZE);
  pthread_mutex_init(&Write_Proto_Mutex, NULL);

  for (uint32 i=0;  i<G.Num_PThreads;  i++) {
    Initialize_Work_Area(thread_wa+i, i);
    thread_wa[i].gkpStore = gkpStore;
  }

  //  Command line options are Lo_Hash_Frag and Hi_Hash_Frag
  //  Command line options are Lo_Old_Frag and Hi_Old_Frag

  if (G.bgnHashID < 1)
    G.bgnHashID = 1;

  if (gkpStore->gkStore_getNumReads() < G.endHashID)
    G.endHashID = gkpStore->gkStore_getNumReads();


  //  Note distinction between the local bgn/end and the global G.bgn/G.end.

  uint32  bgnHashID = G.bgnHashID;
  uint32  endHashID = G.bgnHashID + G.Max_Hash_Strings - 1;  //  Inclusive!

  //  Iterate over read blocks, build a hash table, then search in threads.

  while (bgnHashID < G.endHashID) {
    if (endHashID > G.endHashID)
      endHashID = G.endHashID;

    assert(0          <  bgnHashID);
    assert(bgnHashID  <= endHashID);
    assert(endHashID  <= gkpStore->gkStore_getNumReads());

    //  Load as much as we can.  If we load less than expected, the endHashID is updated to reflect
    //  the last read loaded.

    endHashID = Build_Hash_Index(gkpStore, bgnHashID, endHashID);

    fprintf(stderr, "Index built.\n");

    //  Decide the range of reads to process.  No more than what is loaded in the table.

    if (G.bgnRefID < 1)
      G.bgnRefID = 1;

    if (G.endRefID > gkpStore->gkStore_getNumReads())
      G.endRefID = gkpStore->gkStore_getNumReads();

    G.curRefID = G.bgnRefID;

    //  The old version used to further divide the ref range into blocks of at most
    //  Max_Reads_Per_Batch so that those reads could be loaded into core.  We don't need to do that
    //  anymore.

    G.perThread = (G.endRefID - G.bgnRefID) / G.Num_PThreads / 8;

    if (G.perThread < 1)
      G.perThread = 1;

    fprintf(stderr, "\n");
    fprintf(stderr, "Range: %u-%u.  Store has %u reads.\n",
            G.bgnRefID, G.endRefID, gkpStore->gkStore_getNumReads());
    fprintf(stderr, "Chunk: "F_U32" reads/thread -- (G.endRefID="F_U32" - G.bgnRefID="F_U32") / G.Num_PThreads="F_U32" / 8\n",
            G.perThread, G.endRefID, G.bgnRefID, G.Num_PThreads);

    fprintf(stderr, "\n");
    fprintf(stderr, "Starting "F_U32"-"F_U32" with "F_U32" per thread\n", G.bgnRefID, G.endRefID, G.perThread);
    fprintf(stderr, "\n");

    for (uint32 i=0; i<G.Num_PThreads; i++) {

      //  Initialize each thread, reset the current position.

      thread_wa[i].bgnID = G.curRefID;
      thread_wa[i].endID = thread_wa[i].bgnID + G.perThread - 1;

      G.curRefID = thread_wa[i].endID + 1;

      if (G.endRefID > G.endRefID)
        G.endRefID = G.endRefID;

      int status = pthread_create(thread_id+i, &attr, Process_Overlaps, thread_wa+i);

      if (status != 0)
        fprintf(stderr, "pthread_create error:  %s\n", strerror(status)), exit(1);
    }

    //  The master thread just sits here and waits.

    for (uint32 i=0; i<G.Num_PThreads; i++) {
      int status = pthread_join(thread_id[i], NULL);
      if (status != 0)
        fprintf(stderr, "pthread_join error: %s\n", strerror(status)), exit(1);
    }

    //  Clear out the hash table.  This stuff is allocated in Build_Hash_Index

    delete [] basesData;  basesData = NULL;
    delete [] qualsData;  qualsData = NULL;
    delete [] nextRef;    nextRef   = NULL;

    //  This one could be left allocated, except for the last iteration.

    delete [] Extra_Ref_Space;  Extra_Ref_Space = NULL;  Max_Extra_Ref_Space = 0;

    //  Prepare for another hash table iteration.
    bgnHashID = endHashID + 1;
    endHashID = bgnHashID + G.Max_Hash_Strings - 1;  //  Inclusive!
  }

  pthread_mutex_destroy(&Write_Proto_Mutex);
  pthread_attr_destroy(&attr);

  delete gkpStore;

  for (uint32 i=0;  i<G.Num_PThreads;  i++)
    Delete_Work_Area(thread_wa + i);

  delete [] thread_wa;
  delete [] thread_id;

  return  0;
}





int
main(int argc, char **argv) {
  int  illegal;
  char  * p;

  argc = AS_configure(argc, argv);

  G.initialize();

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      G.Doing_Partial_Overlaps = TRUE;

    } else if (strcmp(argv[arg], "-h") == 0) {
      AS_UTL_decodeRange(argv[++arg], G.bgnHashID, G.endHashID);

    } else if (strcmp(argv[arg], "-H") == 0) {
      AS_UTL_decodeRange(argv[++arg], G.minLibToHash, G.maxLibToHash);

    } else if (strcmp(argv[arg], "-r") == 0) {
      AS_UTL_decodeRange(argv[++arg], G.bgnRefID, G.endRefID);

    } else if (strcmp(argv[arg], "-R") == 0) {
      AS_UTL_decodeRange(argv[++arg], G.minLibToRef, G.maxLibToRef);

    } else if (strcmp(argv[arg], "-k") == 0) {
      arg++;
      if ((isdigit(argv[arg][0]) && (argv[arg][1] == 0)) ||
          (isdigit(argv[arg][0]) && isdigit(argv[arg][1]) && (argv[arg][2] == 0))) {
        G.Kmer_Len = strtoull(argv[arg], NULL, 10);
      } else {
        errno = 0;
        G.Kmer_Skip_File = fopen(argv[arg], "r");
        if (errno)
          fprintf(stderr, "ERROR: Failed to open -k '%s': %s\n", argv[arg], strerror(errno)), exit(1);
      }

    } else if (strcmp(argv[arg], "-l") == 0) {
      G.Frag_Olap_Limit = strtol(argv[++arg], NULL, 10);
      if  (G.Frag_Olap_Limit < 1)
        G.Frag_Olap_Limit = UINT64_MAX;

    } else if (strcmp(argv[arg], "-m") == 0) {
      G.Unique_Olap_Per_Pair = FALSE;
    } else if (strcmp(argv[arg], "-u") == 0) {
      G.Unique_Olap_Per_Pair = TRUE;

    } else if (strcmp(argv[arg], "--hashbits") == 0) {
      G.Hash_Mask_Bits = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--hashstrings") == 0) {
      G.Max_Hash_Strings = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--hashdatalen") == 0) {
      G.Max_Hash_Data_Len = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--hashload") == 0) {
      G.Max_Hash_Load = atof(argv[++arg]);

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
      G.Max_Reads_Per_Batch = strtoul(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--readsperthread") == 0) {
      G.Max_Reads_Per_Thread = strtoul(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-o") == 0) {
      G.Outfile_Name = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      G.Num_PThreads = strtoull(argv[++arg], NULL, 10);


    } else if (strcmp(argv[arg], "-v") == 0) {
      G.Min_Olap_Len = (int) strtol (argv[++arg], & p, 10);

    } else if (strcmp(argv[arg], "-w") == 0) {
      G.Use_Window_Filter = TRUE;

    } else if (strcmp(argv[arg], "-x") == 0) {
      G.Ignore_Clear_Range = TRUE;

    } else if (strcmp(argv[arg], "-z") == 0) {
      G.Use_Hopeless_Check = FALSE;

    } else {
      if (G.Frag_Store_Path == NULL) {
        G.Frag_Store_Path = argv[arg];
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
    if (G.Use_Window_Filter)
      fprintf(stderr, "High error rates requested -- window-filter turned off despite -w flag!\n");
    G.Use_Window_Filter  = FALSE;
    G.Use_Hopeless_Check = FALSE;
  }

  if (G.Max_Hash_Strings == 0)
    fprintf(stderr, "* No memory model supplied; -M needed!\n"), err++;

  if (G.Kmer_Len == 0)
    fprintf(stderr, "* No kmer length supplied; -k needed!\n"), err++;

  if (G.Max_Hash_Strings > MAX_STRING_NUM)
    fprintf(stderr, "Too many strings (--hashstrings), must be less than "F_U64"\n", MAX_STRING_NUM), err++;

  if (G.Outfile_Name == NULL)
    fprintf (stderr, "ERROR:  No output file name specified\n"), err++;

  if ((err) || (G.Frag_Store_Path == NULL)) {
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
    fprintf(stderr, "                     maxreadlen 2048->hashstrings  524288 (default)\n");
    fprintf(stderr, "                     maxreadlen  512->hashstrings 2097152\n");
    fprintf(stderr, "                     maxreadlen  128->hashstrings 8388608\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--readsperbatch n  Force batch size to n.\n");
    fprintf(stderr, "--readsperthread n Force each thread to process n reads.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  Out_BOF = new ovFile(G.Outfile_Name, ovFileFullWrite);

  //  Adjust the number of reads to load into memory at once (for processing, not the hash table),

  if (G.Max_Reads_Per_Batch == 0)
    G.Max_Reads_Per_Batch = (G.Max_Hash_Strings < 100000) ? G.Max_Hash_Strings : 100000;

  //if (Max_Hash_Strings < Max_Reads_Per_Batch)
  //  Max_Reads_Per_Batch = Max_Hash_Strings;

  //  Adjust the number of reads processed per thread.  Default to having four blocks per thread,
  //  but make sure that (a) all threads have work to do, and (b) batches are not minuscule.

  if (G.Max_Reads_Per_Thread == 0)
    G.Max_Reads_Per_Thread = G.Max_Reads_Per_Batch / (4 * G.Num_PThreads);

  if (G.Max_Reads_Per_Thread * G.Num_PThreads > G.Max_Reads_Per_Batch)
    G.Max_Reads_Per_Thread = G.Max_Reads_Per_Batch / G.Num_PThreads + 1;

  if (G.Max_Reads_Per_Thread < 10)
    G.Max_Reads_Per_Thread = 10;

  //  We know enough now to set the hash function variables, and some other random variables.

  HSF1 = G.Kmer_Len - (G.Hash_Mask_Bits / 2);
  HSF2 = 2 * G.Kmer_Len - G.Hash_Mask_Bits;
  SV1  = HSF1 + 2;
  SV2  = (HSF1 + HSF2) / 2;
  SV3  = HSF2 - 2;


  fprintf(stderr, "\n");
  fprintf(stderr, "STRING_NUM_BITS       "F_U32"\n", STRING_NUM_BITS);
  fprintf(stderr, "OFFSET_BITS           "F_U32"\n", OFFSET_BITS);
  fprintf(stderr, "STRING_NUM_MASK       "F_U64"\n", STRING_NUM_MASK);
  fprintf(stderr, "OFFSET_MASK           "F_U64"\n", OFFSET_MASK);
  fprintf(stderr, "MAX_STRING_NUM        "F_U64"\n", MAX_STRING_NUM);
  fprintf(stderr, "\n");
  fprintf(stderr, "Hash_Mask_Bits        "F_U32"\n", G.Hash_Mask_Bits);
  fprintf(stderr, "Max_Hash_Strings      "F_U32"\n", G.Max_Hash_Strings);
  fprintf(stderr, "Max_Hash_Data_Len     "F_U64"\n", G.Max_Hash_Data_Len);
  fprintf(stderr, "Max_Hash_Load         %f\n", G.Max_Hash_Load);
  fprintf(stderr, "Kmer Length           %d\n", G.Kmer_Len);
  fprintf(stderr, "Min Overlap Length    %d\n", G.Min_Olap_Len);
  fprintf(stderr, "MAX_ERRORS            %d\n", MAX_ERRORS);
  fprintf(stderr, "ERRORS_FOR_FREE       %d\n", ERRORS_FOR_FREE);
  fprintf(stderr, "\n");
  fprintf(stderr, "Num_PThreads          "F_U32"\n", G.Num_PThreads);
  fprintf(stderr, "Max_Reads_Per_Batch   "F_U32"\n", G.Max_Reads_Per_Batch);
  fprintf(stderr, "Max_Reads_Per_Thread  "F_U32"\n", G.Max_Reads_Per_Thread);

  assert (8 * sizeof (uint64) > 2 * G.Kmer_Len);


  Bit_Equivalent['a'] = Bit_Equivalent['A'] = 0;
  Bit_Equivalent['c'] = Bit_Equivalent['C'] = 1;
  Bit_Equivalent['g'] = Bit_Equivalent['G'] = 2;
  Bit_Equivalent['t'] = Bit_Equivalent['T'] = 3;

  for  (int i = 0;  i < 256;  i ++) {
    char  ch = tolower ((char) i);

    if  (ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
      Char_Is_Bad[i] = 0;
    else
      Char_Is_Bad[i] = 1;
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "HASH_TABLE_SIZE         "F_U32"\n",     HASH_TABLE_SIZE);
  fprintf(stderr, "sizeof(Hash_Bucket_t)   "F_SIZE_T"\n",  sizeof(Hash_Bucket_t));
  fprintf(stderr, "hash table size:        "F_SIZE_T" MB\n",  (HASH_TABLE_SIZE * sizeof(Hash_Bucket_t)) >> 20);
  fprintf(stderr, "\n");

  Hash_Table       = new Hash_Bucket_t [HASH_TABLE_SIZE];

  fprintf(stderr, "check  "F_SIZE_T" MB\n", (HASH_TABLE_SIZE    * sizeof (Check_Vector_t) >> 20));
  fprintf(stderr, "info   "F_SIZE_T" MB\n", (G.Max_Hash_Strings * sizeof (Hash_Frag_Info_t) >> 20));
  fprintf(stderr, "start  "F_SIZE_T" MB\n", (G.Max_Hash_Strings * sizeof (int64) >> 20));
  fprintf(stderr, "\n");

  Hash_Check_Array = new Check_Vector_t [HASH_TABLE_SIZE];
  String_Info      = new Hash_Frag_Info_t [G.Max_Hash_Strings];
  String_Start     = new int64 [G.Max_Hash_Strings];

  String_Start_Size = G.Max_Hash_Strings;

  memset(Hash_Check_Array, 0, sizeof(Check_Vector_t)   * HASH_TABLE_SIZE);
  memset(String_Info,      0, sizeof(Hash_Frag_Info_t) * G.Max_Hash_Strings);
  memset(String_Start,     0, sizeof(int64)            * G.Max_Hash_Strings);



  OverlapDriver();



  delete [] basesData;
  delete [] qualsData;
  delete [] nextRef;

  delete [] String_Start;
  delete [] String_Info;
  delete [] Hash_Check_Array;
  delete [] Hash_Table;

  fprintf (stderr, " Kmer hits without olaps = "F_S64"\n", Kmer_Hits_Without_Olap_Ct);
  fprintf (stderr, "    Kmer hits with olaps = "F_S64"\n", Kmer_Hits_With_Olap_Ct);
  fprintf (stderr, "  Multiple overlaps/pair = "F_S64"\n", Multi_Overlap_Ct);
  fprintf (stderr, " Total overlaps produced = "F_S64"\n", Total_Overlaps);
  fprintf (stderr, "      Contained overlaps = "F_S64"\n", Contained_Overlap_Ct);
  fprintf (stderr, "       Dovetail overlaps = "F_S64"\n", Dovetail_Overlap_Ct);
  fprintf (stderr, "Rejected by short window = "F_S64"\n", Bad_Short_Window_Ct);
  fprintf (stderr, " Rejected by long window = "F_S64"\n", Bad_Long_Window_Ct);

  delete Out_BOF;

  return(0);
}
