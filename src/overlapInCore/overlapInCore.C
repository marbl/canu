
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

#include "system.H"
#include "strings.H"

#include "overlapInCore.H"

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
uint64  Kmer_Hits_Skipped_Ct = 0;
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

ovFile  *Out_BOF = NULL;



//  Allocate memory for  (* WA)  and set initial values.
//  Set  thread_id  field to  id .
void
Initialize_Work_Area(Work_Area_t *WA, int id, sqStore *readStore, sqCache *readCache) {
  //uint64  allocated = 0;

  WA->String_Olap_Size  = INIT_STRING_OLAP_SIZE;
  WA->String_Olap_Space = new String_Olap_t [WA->String_Olap_Size];

  WA->Match_Node_Size  = INIT_MATCH_NODE_SIZE;
  WA->Match_Node_Space = new Match_Node_t [WA->Match_Node_Size];

  //allocated += WA->String_Olap_Size * sizeof (String_Olap_t);
  //allocated += WA->Match_Node_Size  * sizeof (Match_Node_t);

  WA->status     = 0;
  WA->thread_id  = id;

  WA->readStore = readStore;
  WA->readCache = readCache;

  WA->overlapsLen = 0;
  WA->overlapsMax = 1024 * 1024 / sizeof(ovOverlap);
  WA->overlaps    = new ovOverlap [WA->overlapsMax];

  //allocated += sizeof(ovOverlap) * WA->overlapsMax;

  WA->editDist = new prefixEditDistance(G.Doing_Partial_Overlaps,
                                        G.maxErate,
                                        G.maxErate * G.alignNoise);

  WA->q_diff = new char [AS_MAX_READLEN];
  WA->distinct_olap = new Olap_Info_t [MAX_DISTINCT_OLAPS];
}


void
Delete_Work_Area(Work_Area_t *WA) {
  delete    WA->editDist;
  delete [] WA->String_Olap_Space;
  delete [] WA->Match_Node_Space;
  delete [] WA->overlaps;

  delete [] WA->distinct_olap;
  delete [] WA->q_diff;
}




int
OverlapDriver(void) {

  Work_Area_t    *thread_wa = new Work_Area_t [G.Num_PThreads];

  sqStore        *readStore = new sqStore(G.Frag_Store_Path);
  sqCache        *readCache = new sqCache(readStore);

  Out_BOF = new ovFile(readStore, G.Outfile_Name, ovFileFullWrite);

  fprintf(stderr, "Initializing %u work areas.\n", G.Num_PThreads);

#pragma omp parallel for
  for (uint32 i=0;  i<G.Num_PThreads;  i++)
    Initialize_Work_Area(thread_wa+i, i, readStore, readCache);

  //  Make sure both the hash and reference ranges are valid.

  if (G.bgnHashID < 1)
    G.bgnHashID = 1;

  if (readStore->sqStore_lastReadID() < G.endHashID)
    G.endHashID = readStore->sqStore_lastReadID();

  if (G.bgnRefID < 1)
    G.bgnRefID = 1;

  if (G.endRefID > readStore->sqStore_lastReadID())
    G.endRefID = readStore->sqStore_lastReadID();

  //  Load the reference range into the cache

  fprintf(stderr, "Loading reference reads %u-%u inclusive.\n", G.bgnRefID, G.endRefID);

  readCache->sqCache_loadReads(G.bgnRefID, G.endRefID, true);

  //  Note distinction between the local bgn/end and the global G.bgn/G.end.

  uint32  bgnHashID = G.bgnHashID;
  uint32  endHashID = G.endHashID;

  //  Iterate over read blocks, build a hash table, then search in threads.

  while (bgnHashID < G.endHashID) {
    if (endHashID > G.endHashID)
      endHashID = G.endHashID;

    assert(0          <  bgnHashID);
    assert(bgnHashID  <= endHashID);
    assert(endHashID  <= readStore->sqStore_lastReadID());

    //  Load as much as we can.  If we load less than expected, the endHashID is updated to reflect
    //  the last read loaded.

    endHashID = Build_Hash_Index(readStore, bgnHashID, endHashID);

    //  Decide the range of reads to process.  No more than what is loaded in the table.

    G.curRefID = G.bgnRefID;

    //  The old version used to further divide the ref range into blocks of at most
    //  Max_Reads_Per_Batch so that those reads could be loaded into core.  We don't
    //  need to do that anymore.

    G.perThread = 1 + (G.endRefID - G.bgnRefID) / G.Num_PThreads / 8;

    fprintf(stderr, "\n");
    fprintf(stderr, "Range: %u-%u.  Store has %u reads.\n",
            G.bgnRefID, G.endRefID, readStore->sqStore_lastReadID());
    fprintf(stderr, "Chunk: " F_U32 " reads/thread -- (G.endRefID=" F_U32 " - G.bgnRefID=" F_U32 ") / G.Num_PThreads=" F_U32 " / 8\n",
            G.perThread, G.endRefID, G.bgnRefID, G.Num_PThreads);

    fprintf(stderr, "\n");
    fprintf(stderr, "Starting " F_U32 "-" F_U32 " with " F_U32 " per thread\n", G.bgnRefID, G.endRefID, G.perThread);
    fprintf(stderr, "\n");

    //  Initialize each thread, reset the current position.  curRefID and endRefID are updated, this
    //  cannot be done in the parallel loop!

    for (uint32 i=0; i<G.Num_PThreads; i++) {
      thread_wa[i].bgnID = G.curRefID;
      thread_wa[i].endID = thread_wa[i].bgnID + G.perThread - 1;

      G.curRefID = thread_wa[i].endID + 1;  //  Global value updated!
    }

#pragma omp parallel for
    for (uint32 i=0; i<G.Num_PThreads; i++)
      Process_Overlaps(thread_wa + i);

    //  Clear out the hash table.  This stuff is allocated in Build_Hash_Index

    delete [] basesData;  basesData = NULL;
    delete [] nextRef;    nextRef   = NULL;

    //  This one could be left allocated, except for the last iteration.

    delete [] Extra_Ref_Space;  Extra_Ref_Space = NULL;  Max_Extra_Ref_Space = 0;

    //  Prepare for another hash table iteration.
    bgnHashID = endHashID + 1;
    endHashID = G.endHashID;
  }

  delete Out_BOF;

  delete readCache;

  delete readStore;

  for (uint32 i=0;  i<G.Num_PThreads;  i++)
    Delete_Work_Area(thread_wa + i);

  delete [] thread_wa;

  return  0;
}





int
main(int argc, char **argv) {
  int  illegal;

  argc = AS_configure(argc, argv);

  G.initialize();  //  Probably redundant with the call in the constructor, but doesn't hurt.

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-partial") == 0) {
      G.Doing_Partial_Overlaps = true;

    } else if (strcmp(argv[arg], "-h") == 0) {
      decodeRange(argv[++arg], G.bgnHashID, G.endHashID);

    } else if (strcmp(argv[arg], "-H") == 0) {
      decodeRange(argv[++arg], G.minLibToHash, G.maxLibToHash);

    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], G.bgnRefID, G.endRefID);

    } else if (strcmp(argv[arg], "-R") == 0) {
      decodeRange(argv[++arg], G.minLibToRef, G.maxLibToRef);

    } else if (strcmp(argv[arg], "-k") == 0) {
      arg++;

      if (fileExists(argv[arg]) == true)
        G.kmerSkipFileName = argv[arg];
      else
        G.Kmer_Len = strtouint32(argv[arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      G.Frag_Olap_Limit = strtol(argv[++arg], NULL, 10);
      if  (G.Frag_Olap_Limit < 1)
        G.Frag_Olap_Limit = UINT64_MAX;

    } else if (strcmp(argv[arg], "-m") == 0) {
      G.Unique_Olap_Per_Pair = false;
    } else if (strcmp(argv[arg], "-u") == 0) {
      G.Unique_Olap_Per_Pair = true;

    } else if (strcmp(argv[arg], "--hashbits") == 0) {
      G.Hash_Mask_Bits = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--hashdatalen") == 0) {
      G.Max_Hash_Data_Len = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "--hashload") == 0) {
      G.Max_Hash_Load = atof(argv[++arg]);

#if 0
    //  This should still work, but not useful unless String_Ref_t is
    //  changed to uint32.
    //
    //fprintf(stderr, "--maxreadlen n     For batches with all short reads, pack bits differently to\n");
    //fprintf(stderr, "                   process more reads per batch.\n");
    //fprintf(stderr, "                     all reads must be shorter than n\n");
    //fprintf(stderr, "                     limited to 2^(30-m) reads\n");
    //fprintf(stderr, "                   Common values:\n");
    //fprintf(stderr, "                     maxreadlen 2048->hashstrings  524288 (default)\n");
    //fprintf(stderr, "                     maxreadlen  512->hashstrings 2097152\n");
    //fprintf(stderr, "                     maxreadlen  128->hashstrings 8388608\n");
    //fprintf(stderr, "\n");

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
#endif

    } else if (strcmp(argv[arg], "-o") == 0) {
      G.Outfile_Name = argv[++arg];

    } else if (strcmp(argv[arg], "-s") == 0) {
      G.Outstat_Name = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      G.Num_PThreads = setNumThreads(argv[++arg]);


    } else if (strcmp(argv[arg], "--minlength") == 0) {
      G.Min_Olap_Len = strtol (argv[++arg], NULL, 10);
    } else if (strcmp(argv[arg], "--minkmers") == 0) {
      G.Filter_By_Kmer_Count = 1;
    } else if (strcmp(argv[arg], "--maxerate") == 0) {
      G.maxErate = strtof(argv[++arg], NULL);
    } else if (strcmp(argv[arg], "--alignnoise") == 0) {
      G.alignNoise = strtof(argv[++arg], NULL);

    } else if (strcmp(argv[arg], "-z") == 0) {
      G.Use_Hopeless_Check = false;

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
  if (G.maxErate > 0.06)
    G.Use_Hopeless_Check = false;

  if (G.Kmer_Len == 0)
    fprintf(stderr, "* No kmer length supplied; -k needed!\n"), err++;

  if (G.Outfile_Name == NULL)
    fprintf (stderr, "ERROR:  No output file name specified\n"), err++;

   // if we were asked to use the k-mer filter, initialize its value now that we know the parameters
   if (G.Filter_By_Kmer_Count == 1)
      G.Filter_By_Kmer_Count = int(floor(exp(-1.0 * (double)G.Kmer_Len * G.maxErate) * (G.Min_Olap_Len - G.Kmer_Len + 1)));

  if ((err) || (G.Frag_Store_Path == NULL)) {
    fprintf(stderr, "USAGE:  %s [options] <seqStorePath>\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "-b <fn>     in contig mode, specify the output file\n");
    fprintf(stderr, "-c          contig mode.  Use 2 frag stores.  First is\n");
    fprintf(stderr, "            for reads; second is for contigs\n");
    fprintf(stderr, "-partial    do partial overlaps\n");
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
    fprintf(stderr, "-P          write protoIO output (if not -partial)\n");
    fprintf(stderr, "-r <range>  specify old fragments to overlap\n");
    fprintf(stderr, "-t <n>      use <n> parallel threads\n");
    fprintf(stderr, "-u          allow only 1 overlap per oriented fragment pair\n");
    fprintf(stderr, "-z          skip the hopeless check (also skipped at > 0.06)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--maxerate <n>     only output overlaps with fraction <n> or less error (e.g., 0.06 == 6%%)\n");
    fprintf(stderr, "--minlength <n>    only output overlaps of <n> or more bases\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--hashbits n       Use n bits for the hash mask.\n");
    fprintf(stderr, "--hashdatalen n    Load at most n bytes into the hash table at one time.\n");
    fprintf(stderr, "--hashload f       Load to at most 0.0 < f < 1.0 capacity (default 0.7).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--readsperbatch n  Force batch size to n.\n");
    fprintf(stderr, "--readsperthread n Force each thread to process n reads.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  //  We know enough now to set the hash function variables, and some other random variables.
    
  G.Max_Hash_Data_Len += AS_MAX_READLEN;
  HSF1 = G.Kmer_Len - (G.Hash_Mask_Bits / 2);
  HSF2 = 2 * G.Kmer_Len - G.Hash_Mask_Bits;
  SV1  = HSF1 + 2;
  SV2  = (HSF1 + HSF2) / 2;
  SV3  = HSF2 - 2;

  //  Log parameters.

#if 0
  fprintf(stderr, "\n");
  fprintf(stderr, "STRING_NUM_BITS          " F_U32 "\n", STRING_NUM_BITS);
  fprintf(stderr, "OFFSET_BITS              " F_U32 "\n", OFFSET_BITS);
  fprintf(stderr, "STRING_NUM_MASK          " F_U64 "\n", STRING_NUM_MASK);
  fprintf(stderr, "OFFSET_MASK              " F_U64 "\n", OFFSET_MASK);
  fprintf(stderr, "MAX_STRING_NUM           " F_U64 "\n", MAX_STRING_NUM);
  fprintf(stderr, "\n");
  fprintf(stderr, "Hash_Mask_Bits           " F_U32 "\n", G.Hash_Mask_Bits);
  fprintf(stderr, "\n");
  fprintf(stderr, "bgnHashID                " F_U32 "\n", G.bgnHashID);
  fprintf(stderr, "bgnHashID                " F_U32 "\n", G.endHashID);
  fprintf(stderr, "\n");
  fprintf(stderr, "Max_Hash_Data_Len        " F_U64 "\n", G.Max_Hash_Data_Len);
  fprintf(stderr, "Max_Hash_Load            %f\n", G.Max_Hash_Load);
  fprintf(stderr, "Kmer Length              " F_U64 "\n", G.Kmer_Len);
  fprintf(stderr, "Min Overlap Length       %d\n", G.Min_Olap_Len);
  fprintf(stderr, "Max Error Rate           %f\n", G.maxErate);
  fprintf(stderr, "Min Kmer Matches         " F_U64 "\n", G.Filter_By_Kmer_Count);
  fprintf(stderr, "\n");
  fprintf(stderr, "Num_PThreads             " F_U32 "\n", G.Num_PThreads);
  fprintf(stderr, "\n");
  fprintf(stderr, "sizeof(Hash_Bucket_t)    " F_U64 "\n",     (uint64)sizeof(Hash_Bucket_t));
  fprintf(stderr, "sizeof(Check_Vector_t)   " F_U64 "\n",     (uint64)sizeof(Check_Vector_t));
  fprintf(stderr, "sizeof(Hash_Frag_Info_t) " F_U64 "\n",     (uint64)sizeof(Hash_Frag_Info_t));
  fprintf(stderr, "\n");
  fprintf(stderr, "HASH_TABLE_SIZE          " F_U64 "\n",     HASH_TABLE_SIZE);
  fprintf(stderr, "\n");
  fprintf(stderr, "hash table size:         " F_U64    " MB\n", (HASH_TABLE_SIZE * sizeof(Hash_Bucket_t)) >> 20);
  fprintf(stderr, "hash check array         " F_U64    " MB\n", (HASH_TABLE_SIZE    * sizeof (Check_Vector_t))   >> 20);
  fprintf(stderr, "string info              " F_SIZE_T " MB\n", ((G.endHashID - G.bgnHashID + 1) * sizeof (Hash_Frag_Info_t)) >> 20);
  fprintf(stderr, "string start             " F_SIZE_T " MB\n", ((G.endHashID - G.bgnHashID + 1) * sizeof (int64))            >> 20);
  fprintf(stderr, "\n");
#endif

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

  Hash_Table       = new Hash_Bucket_t    [HASH_TABLE_SIZE];
  Hash_Check_Array = new Check_Vector_t   [HASH_TABLE_SIZE];
  String_Info      = new Hash_Frag_Info_t [G.endHashID - G.bgnHashID + 1];
  String_Start     = new int64            [G.endHashID - G.bgnHashID + 1];

  String_Start_Size = G.endHashID - G.bgnHashID + 1;

  memset(Hash_Check_Array, 0, sizeof(Check_Vector_t)   * HASH_TABLE_SIZE);
  memset(String_Info,      0, sizeof(Hash_Frag_Info_t) * (G.endHashID - G.bgnHashID + 1));
  memset(String_Start,     0, sizeof(int64)            * (G.endHashID - G.bgnHashID + 1));



  OverlapDriver();



  delete [] basesData;
  delete [] nextRef;

  delete [] String_Start;
  delete [] String_Info;
  delete [] Hash_Check_Array;
  delete [] Hash_Table;

  FILE *stats = stderr;

  if (G.Outstat_Name != NULL) {
    errno = 0;
    stats = fopen(G.Outstat_Name, "w");
    if (errno) {
      fprintf(stderr, "WARNING: failed to open '%s' for writing: %s\n", G.Outstat_Name, strerror(errno));
      stats = stderr;
    }
  }

  fprintf(stats, " Kmer hits without olaps = " F_S64 "\n", Kmer_Hits_Without_Olap_Ct);
  fprintf(stats, "    Kmer hits with olaps = " F_S64 "\n", Kmer_Hits_With_Olap_Ct);
  //fprintf(stats, "      Kmer hits below %u = " F_S64 "\n", G.Filter_By_Kmer_Count, Kmer_Hits_Skipped_Ct);
  fprintf(stats, "  Multiple overlaps/pair = " F_S64 "\n", Multi_Overlap_Ct);
  fprintf(stats, " Total overlaps produced = " F_S64 "\n", Total_Overlaps);
  fprintf(stats, "      Contained overlaps = " F_S64 "\n", Contained_Overlap_Ct);
  fprintf(stats, "       Dovetail overlaps = " F_S64 "\n", Dovetail_Overlap_Ct);
  fprintf(stats, "Rejected by short window = " F_S64 "\n", Bad_Short_Window_Ct);
  fprintf(stats, " Rejected by long window = " F_S64 "\n", Bad_Long_Window_Ct);

  merylutil::closeFile(stats, G.Outstat_Name);

  fprintf(stderr, "Bye.\n");

  return(0);
}
