//  Don't define this (unless your name starts with L or B).
//
//  Changes to here, sim4db.H and exon_cores.H
//
//#define INTERSPECIES

#define  DIST_CUTOFF       3

#define  MIN_INTRON        30
#define  SHORT_INTRON      50
#define  LONG_INTRON       20000
#define  SHORT_EXON        40
#define  MAX_GRINIT        500
#define  MAX_SLIDE         15
#define  MAX_INTERNAL_GAP  50

#define  DEFAULT_DRANGE    10
#define  DEFAULT_WEIGHT    100
#define  DEFAULT_RELINK_WEIGHT    500

#define  DEFAULT_K         16
#define  DEFAULT_C         12

#define  min(x,y)        ((x>y) ? (y):(x))
#define  max(x,y)        ((x<y) ? (y):(x))

#define  START_SIG       ((G_score >= abs(C_score)) ? "GT" : "CT") 
#define  END_SIG         ((G_score >= abs(C_score)) ? "AG" : "AC")

#define  DELETE          1
#define  INSERT          2
#define  SUBSTITUTE      3
#define  INTRON          4
#define  O_INTRON        5

#define  SPLICE_ORIGINAL         0
#define  SPLICE_GENESPLICER      1
#define  SPLICE_GLIMMER          2
#define  DEFAULT_SPLICE_MODEL    SPLICE_ORIGINAL

#define  DEFAULT_SPACED_SEED        "1xx1011011011xx11"
#define  DEFAULT_SPACED_SEED_INT    "10011010100011"
#define  DEFAULT_SPACED_SEED_EXT    "1101100011010111"
#define  SEED_SPAN               40   /* 22 */

