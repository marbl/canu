#define  BAD_WINDOW_LEN           50
//  Length of window in which to look for clustered errors
//  to invalidate an overlap

#define  BAD_WINDOW_VALUE         (8 * QUALITY_CUTOFF)
//  This many or more errors in a window of  BAD_WINDOW_LEN
//  invalidates an overlap

#define  DEFAULT_HI_HIT_LIMIT    INT_MAX
//  Any kmer in hash table with at least this many hits
//  cannot initiate an overlap.  Can be changed on the command
//  line with  -K  option

#define  DELETED_FRAG            2
//  Indicates fragment was marked as deleted in the store

#define  DISPLAY_WIDTH           60
//  Number of characters per line when displaying sequences

#define  ENTRIES_PER_BUCKET      21
//  In main hash table.  Recommended values are 21, 31 or 42
//  depending on cache line size.


#define  HASH_EXPANSION_FACTOR   1.4
//  Hash table size is >= this times  MAX_HASH_STRINGS

#define  HASH_TABLE_SIZE         (1 + HASH_MASK)
//  Number of buckets in hash table

#define  HOPELESS_MATCH          90
//  A string this long or longer without an exact kmer
//  match is assumed to be hopeless to find a match
//  within the error threshold

#define  IID_GAP_LIMIT           100
//  When using a list of fragment IID's, gaps between
//  IID's this long or longer force a new load partial
//  store to be done

#define  INIT_SCREEN_MATCHES     50
//  Initial number of screen-match entries per fragment

#define  K_MER_STEP          1
//  1 = every k-mer in search
//  2 = every other k-mer
//  3 = every third, ...
//  Used to skip some k-mers in finding matches

#define  MAX_BRANCH_COUNT        UCHAR_MAX
//  The largest branch point count before an overflow

#define  MAX_DISTINCT_OLAPS      3
//  Most possible really different overlaps (i.e., not
//  just shifts from periodic regions) between 2 fragments
//  in a given orientation.  For fragments of approximately
//  same size, should never be more than 2.

#define  MAX_NAME_LEN            500
//  Longest file name allowed

#define  MIN_CALC_KMER           4
//  When calculating the  Hi_Hit_Limit  based on genome length, etc,
//  don't set it below this

#define  MIN_INTERSECTION        10
//  Minimum length of match region to be worth reporting

#define  MIN_OLAP_OUTSIDE_SCREEN 30
//  Minimum number of bases outside of screened regions
//  to be a reportable overlap.  Entire overlap (including screened
//  portion) must still be  MIN_OLAP_LEN .

#define  OUTPUT_OVERLAP_DELTAS   0
//  If true include delta-encoding of overlap alignment
//  in overlap messages.  Otherwise, omit them.
//  As of 6 Oct 2008, support for overlap deltas has been removed.
//  However, there are enough remnants in AS_MSG to output them.
//  Just enabling OUTPUT_OVERLAP_DELTAS will not compile; see
//  AS_MSG_USE_OVL_DELTA in AS_MSG.

#define  QUALITY_CUTOFF          20
//  Regard quality values higher than this as equal to this
//  for purposes of finding bad windows

#define  SCRIPT_NAME             "lsf-ovl"
//  Default name of script produced by  make-ovl-script

#define  SHIFT_SLACK  1
// Allow to be off by this many bases in combining/comparing alignments

#define  STRING_OLAP_SHIFT       8
//  To compute hash function into the String_Olap hash table.

#define  STRING_OLAP_MODULUS     (1 << STRING_OLAP_SHIFT)
//  The size of the  String_Olap  hash table.  The rest of
//  the space is used for chaining.  This number should be
//  relatively small to reflect the number of fragments a
//  given fragment has exact matches with.

#define  STRING_OLAP_MASK        (STRING_OLAP_MODULUS - 1)
//  To compute hash function into the String_Olap hash table.

#define  THREAD_STACKSIZE        (16 * 512 * 512)
//  The amount of stack space to allocate to each thread.

#define  MAX_EXTRA_SUBCOUNT        (AS_MAX_READLEN / G.Kmer_Len)
