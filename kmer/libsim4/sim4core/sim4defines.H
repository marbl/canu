
//  Don't define this (unless your name starts with L or B).
//
//  Changes to here, sim4db.H and exon_cores.H
//
//#define INTERSPECIES


//#define  DEFAULT_NUM_I     3       
#define  DIST_CUTOFF       3
//#define  DEFAULT_MIN_COV   (.8)

#define  MIN_INTRON        30
//#define  MAX_INTRON        20000
#define  MAX_GRINIT        500
#define  MAX_INTERNAL_GAP  50
//#define  LL                60

#define  DEFAULT_DRANGE    10
#define  DEFAULT_WEIGHT    100
#define  DEFAULT_RELINK_WEIGHT    500
#define  DEFAULT_W         12

#ifdef INTERSPECIES
#define  DEFAULT_X         12
#else
#define  DEFAULT_X         12
#endif

#define  DEFAULT_K         16
#define  DEFAULT_C         12

#ifdef INTERSPECIES
#define  P                 (.25)
#else
#define  P                 (.2)
#endif

#define  min(x,y)        ((x>y) ? (y):(x))
#define  max(x,y)        ((x<y) ? (y):(x))
#define  START_SIG       ((G_score >= abs(C_score)) ? "GT" : "CT") 
#define  END_SIG         ((G_score >= abs(C_score)) ? "AG" : "AC")

#ifdef INTERSPECIES
#define  MATCH           1
#define  MISMATCH       (-3)
#else
#define  MATCH           1
#define  MISMATCH       (-5)
#endif

#define  L               8

#define  DELETE          1
#define  INSERT          2
#define  SUBSTITUTE      3
#define  INTRON          4
#define  O_INTRON        5

#ifdef FALSE
#undef FALSE
#endif 

#ifdef TRUE
#undef TRUE
#endif
