#ifndef SIM4_H
#define SIM4_H

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199506L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <assert.h>
#include <sys/types.h>
#include <ctype.h>
#include <errno.h>

#include "util++.H"

#include "ckalloc.H"
#include "sim4parameters.H"
#include "sim4command.H"
#include "exon.H"
#include "mspManager.H"
#include "sim4defines.H"

#include "../sim4polish/sim4polish.h"
#include "../sim4polish/sim4polishList.H"
#include "../sim4polish/sim4polishBuilder.H"

//
//  A wrapper for the sim4 subroutines
//





//  Define this to get a detailed report on hit extension (in mspManager.C)
//
//#define DEBUG_EXTENSION

//  Define this to test the seeds when doing extension -- it tests
//  if the seed is infact an exact match (in mspManager.C)
//
#define TEST_SEEDS_IN_EXTENSION


//  Debug of sim4b1.C and related.
//
//  SHOW_PROGRESS -- write the progress of Sim4::SIM4 to stderr
//  DEBUG_EXONS   -- dump the exons at various places
//
//#define SHOW_PROGRESS
//#define DEBUG_EXONS


//  Show what is being done for external seeding in sim4string.C
//
//#define SHOW_EXTERNAL_SEEDING

enum { INIT = 0, PERM = 1, TEMP = 2};
enum { FWD = 0, BWD = 1, BOTH = 2 };
enum { OK = 0, FREE_START = 1, FREE_END = 2, FREE_BOTH_ENDS = 3};

#define GEN_LOG4_ENTRIES 45
#define CDNA_LOG4_ENTRIES 25

#define HASH_SIZE 32767        /* 2**15 - 1 */

class Sim4 {
private:

struct coords {
  int pos1;
  int pos2;
};

/* used only in the alignment stage */
struct edit_script {
  int                 op_type;    /* SUB, INS, or DEL */
  int                 num;        /* Number of operations */
  struct edit_script *next;
};

struct edit_script_list {
  int   offset1, offset2;
  int   len1, len2; 
  int    score;
  struct edit_script *script;
  struct edit_script_list *next_script;
};

struct splice_t {
  int               xs;
  int               xe;
  int               ys;
  int               ye;
  int               score;
  int               type;
  struct splice_t  *next;
};

struct sim4_stats_t {
  int    internal;
  int    icoverage;
  int    numberOfMatches;
  int    numberOfNs;
  int    orientation;
  int    percentID;
  bool   tooManyMSPs;
};

struct hash_node {
  int                 ecode;              /* integer encoding of the word */
  int                 pos;                /* positions where word hits query sequence */
  struct hash_node   *link; /* next word with same last 7.5 letters */
};

struct hash_table {
  struct hash_node   *table[HASH_SIZE+1];
  int                *nextPos;
  struct hash_node   *nodes;
  int                 nodesused;
};




private:
  mspManager        _mspManager;

  sim4parameters   *globalParams;

  char       *_genSeq;
  char       *_estSeq;

  int         _genLen;
  int         _estLen;

  int         encoding[256];

  coords      last_GT;
  coords      last_CT;
  coords      last_AG;
  coords      last_AC;

  //sim4_args_t rs;

  int         wordExtensionAllowance;

  int         wordSize;
  int         wordSizeInt;
  int         wordSizeExt;

  int         mspThreshold1;
  int         mspThreshold2;

  int         G_score;
  int         C_score;
  Exon       *exon_list;

  int mask;

  struct hash_table  *hashtable;
  struct hash_table   phashtable;
  struct hash_table   thashtable;

  //  For high-accuracy switch
  //
  bool _accurateSequences;


#ifdef __APPLE__
  bool   _pad[7];
#endif



  void adjustBoundariesOfMarginalExons(Exon    *Lblock);
  void findLastAGandAC(Exon *tmp_block1);
  void findLastGTandCT(Exon *tmp_block);




  //  Functions from Xtend1.c
  //
  void free_coords(coords ***val, int size) { 
    for(int i=0; i<size; ++i) 
      ckfree(val[i]);
    ckfree(val);
  };

  int    EXTEND_BW(char *a, char *b, int c, int d, int e, int f, int *g, int *h) {
    if (_accurateSequences) {
      fprintf(stderr, "extend_bw -- Xextend_bw NOT OPTIMIZED\n");
      return(Xextend_bw(a,b,c,d,e,f,g,h));
    } else {
      return(extend_bw(a,b,c,d,e,f,g,h));
    }
  };

  int    EXTEND_FW(char *a, char *b, int c, int d, int e, int f, int *g, int *h) {
    if (_accurateSequences) {
      fprintf(stderr, "extend_bw -- Xextend_fw NOT OPTIMIZED\n");
      return(Xextend_fw(a,b,c,d,e,f,g,h));
    } else {
      return(extend_fw(a,b,c,d,e,f,g,h));
    }
  };

  int    Xextend_bw(char *,char *,int,int,int,int,int *,int *);
  int    Xextend_fw(char *,char *,int,int,int,int,int *,int *);




  //  Functions from align.c
  //
  void align_path(int,int,int,int,int,edit_script**,edit_script**);
  int  align_get_dist(int, int, int, int, int);
  void Condense_script(edit_script *);
  void Condense_both_Ends(edit_script **, edit_script **, edit_script **);
  void Flip_script(struct edit_script **);
#ifdef AUXUTILS
  void Reverse_script(edit_script *);
  void Print_script(edit_script *head, int M, int N);
#endif
  void S2A(edit_script *, int *);
  void align_reverse(int *);
  void Free_script(edit_script *);

  int snake(int k, int x, int endx, int endy) {
    int y;

    if (x<0)
      return x;
    y = x+k;
    while (x<endx && y<endy && _genSeq[x] == _estSeq[y]) {
      ++x;
      ++y;
    }
    return x;
  };


  int rsnake(int k, int x, int startx, int starty, int m) {
    int y;

    if (x>m)
      return x;
    y = x+k;
    while (x>startx && y>starty && _genSeq[x-1] == _estSeq[y-1]) {
      --x;
      --y;
    }
    return x;
  };


  bool  checkExonsForOverlaps(Exon *theExons);
  void  appendExons(sim4polishBuilder &B, Exon *theExon);
  void  maskExonsFromSeeds(sim4command *cmd, Exon *theExons);
  void  maskExonsFromGenomic(Exon *theExons, char *f, char *r, int l);

  void    IDISPLAY(sim4polishBuilder &builder,
                   char *aString,
                   char *bString,
                   char A[],
                   char B[],
                   int M,
                   int N,
                   int S[],
                   int AP,
                   int BP,
                   int est_strand,
                   Exon *exons);
  void    appendAlignments(sim4polishBuilder &builder,
                           char *seq1,
                           char *seq2,
                           int len1,
                           int len2, 
                           edit_script_list *Aligns,
                           Exon *Exons, 
                           int match_ori);








  //  Functions from sim4b1.h
  //
  int good_ratio(int);
  void flip_list(Exon **,Exon **);
  void free_align(edit_script_list *);
  void complement_exons(Exon **,int,int);


  void bld_table(char *,int,int,int);

  int SIM4_block1(Exon*  &Lblock,
                  Exon*  &tmp_block,
                  Exon*  &tmp_block1);
  int SIM4_block2(Exon*   &tmp_Lblock,
                  Exon*   &tmp_Rblock,
                  Exon*   &tmp_block,
                  Exon*   &tmp_block1);
  int SIM4_block3(bool     good_match,
                  Exon*   &tmp_Lblock,
                  Exon*   &tmp_Rblock,
                  Exon*   &tmp_block,
                  Exon*   &tmp_block1);
  int SIM4_block4(bool     good_match,
                  Exon*   &tmp_Lblock,
                  Exon*   &tmp_Rblock,
                  Exon*   &tmp_block,
                  Exon*   &tmp_block1);

  struct edit_script_list *SIM4(int          *,
                                Exon*        *,
                                int          *pA,
                                int          *pT,
                                sim4_stats_t *);

  void   merge(Exon **,Exon **); 
  bool   get_sync_flag(Exon *, Exon *, int);
  void   slide_intron(int w, Exon **, sim4_stats_t *);
  void   sync_slide_intron(int w, Exon **, sim4_stats_t *);
  void   wobble(Exon **,Exon **,const char *,const char *,char *seq1);
  Exon  *bmatch(char *,char *,int,int,int,int);
  Exon  *fmatch(char *,char *,int,int,int,int);
  void   compact_list(Exon **Lblock, Exon **Rblock);
  int    resolve_overlap(Exon *,Exon *,char *);
  int    greedy(char *,char *,int,int,int,int,Exon **, Exon **);
  int    extend_bw(char *,char *,int,int,int,int,int *,int *);
  int    extend_fw(char *,char *,int,int,int,int,int *,int *);

  void   pluri_align(int *dist_ptr,
                     Exon *lblock,
                     struct edit_script_list **Aligns,
                     sim4_stats_t *st);
  void   updateStatistics(Exon *theExon,
                          sim4_stats_t *st);

  void   get_stats(Exon *,sim4_stats_t *); 
  int    get_edist(int,int,int,int,char *,char *);

  void   add_word(int,int);
  int     extend_hit(int,int,const char *const,const char * const,int,int,int,int);
  void   search(char *,char *,int,int,int); 
  void   exon_cores(char*,char*,int,int,int,int,int,int,int,int);

  Exon *find_previous(Exon *head, Exon *target) {
    while (head && (head->next_exon != target))
      head = head->next_exon; 
    return(head);
  };

  bool   get_match_quality(Exon *,Exon *,sim4_stats_t *,int);

  //  splice.h
  //
  void splice_donor(char *xseq, char *yseq, int M, int N, int *gt_score,
                    int *ct_score, int **max_Gf, int **max_Cf,
                    int **start_Gi, int **start_Ci);
  void splice_donor_uni(char *xseq, char *yseq, int M, int N,
                        int *It_score, int **max_IF, int **end_Ii);
  void splice_acceptor(char *xseq, char *yseq, int M, int N,
                       int *ag_score, int *ac_score, int **max_Gb,  
                       int **max_Cb, int **end_Gi, int **end_Ci);
  void splice_acceptor_uni(char *xseq, char *yseq, int M, int N,
                           int *aI_score, int **max_Ib, int **end_Ii);
  int  stepct(int n) {
    if (n<0) fatal("splice.c: Negative value in stepct().");
    if (n<=4) return((int)1.6*8);
    if (n<=8) return((int)1.6*9);
    if (n<=12) return((int)1.6*10);
    return((int)1.6*11);
  };

  splice_t *new_splice(char,int,int,int,int,int,splice_t *);
  void splice(char *in_seqx, int ls, int us, int le, int ue, 
              char *in_seqy, int ys, int ye, splice_t **g, splice_t **c, int ori);






  //  Functions from misc.h
  //
  void fatal(const char *msg) {
    fflush(stdout);
    fprintf(stderr, "%s\n", msg);
    fflush(stderr);
    kill(getpid(), SIGKILL);
  };


  //  Poly-A/T masking stuff
  //
#define T_ONLY    1
#define A_ONLY    2
#define BOTH_AT   3

  void get_polyAT(char *seq,
                  int    len,
                  int   *pA,
                  int   *pT,
                  int    flag=BOTH_AT);

  //void remove_poly(struct edit_script_list **,Exon *,char *,char *,int,int *,int *);


  void remove_polyA_back(struct edit_script_list **,Exon *,char *,char*,int,int *);
  void remove_polyT_front(struct edit_script_list **,Exon *,char *,char*,int *);
  void trim_polyT_align(struct edit_script_list **,Exon **,const int,int *,char *,char *);
  void trim_polyA_align(struct edit_script_list **,Exon *,Exon **,const int,int *,char *,char *);



/* reverse a list of edit script chains */
void
script_flip_list(edit_script_list **left) {  
   edit_script_list  *ep, *ahead, *behind;
       
   ahead = *left;
   ep = NULL;
   while (ahead!=NULL) {
          behind = ep; ep = ahead;
          ahead = ahead->next_script;
          ep->next_script = behind;
  }  
  *left = ep;
}    



public:
  Sim4(sim4parameters *p) {
    globalParams = p;

    _genSeq = 0L;
    _estSeq = 0L;
    _genLen = 0;
    _estLen = 0;

    for (u32bit i=256; i;)
      encoding[--i] = -1;
    encoding[(int)'A'] = encoding[(int)'a'] = 0;
    encoding[(int)'C'] = encoding[(int)'c'] = 1;
    encoding[(int)'G'] = encoding[(int)'g'] = 2;
    encoding[(int)'T'] = encoding[(int)'t'] = 3;

    last_GT.pos1 = last_GT.pos2 = 0;
    last_CT.pos1 = last_CT.pos2 = 0;
    last_AG.pos1 = last_AG.pos2 = 0;
    last_AC.pos1 = last_AC.pos2 = 0;

    wordExtensionAllowance   = 12;

    _mspManager.setLimits(globalParams->_mspLimitAbsolute, globalParams->_mspLimitPercent);
    _mspManager.setParameters(globalParams->_match,
                              globalParams->_mismatch,
                              globalParams->_percentError,
                              wordExtensionAllowance);

    //  wordSize      -- for finding initial seeds
    //  wordSizeInt -- for extending seeds between seeds
    //  wordSizeExt -- for extending seeds on the ends
    //
    //  If sim4parameters defined wordSizeInt or Ext use that,
    //  otherwise, use the original method.
    //
    wordSize      = globalParams->_wordSize;
    wordSizeInt   = min(8,  globalParams->_wordSize);
    wordSizeExt   = min(10, globalParams->_wordSize);

    if (globalParams->_wordSizeInt)  wordSizeInt = globalParams->_wordSizeInt;
    if (globalParams->_wordSizeExt)  wordSizeExt = globalParams->_wordSizeExt;

    mspThreshold1 = globalParams->_mspThresh1;  //K;
    mspThreshold2 = globalParams->_mspThresh2;  //C;

    G_score = 0;
    C_score = 0;

    exon_list = 0L;

    mask = 0;

    hashtable = 0L;

    phashtable.nextPos   = 0L;
    phashtable.nodes     = 0L;
    phashtable.nodesused = 0;

    thashtable.nextPos   = 0L;
    thashtable.nodes     = 0L;
    thashtable.nodesused = 0;

    _accurateSequences = false;
  };

  ~Sim4() {
    delete [] phashtable.nextPos;
    delete [] phashtable.nodes;

    delete [] thashtable.nextPos;
    delete [] thashtable.nodes;

#ifdef WITH_PALLOC
    pfree();
#endif

  };

  sim4polishList *run(sim4command *cmd);
};


#endif  //  SIM4_H
