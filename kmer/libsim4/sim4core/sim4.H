#ifndef SIM4_H
#define SIM4_H

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199506L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
#include <unistd.h>
#include <signal.h>
#include <assert.h>
#include <sys/types.h>
#include <ctype.h>
#include <errno.h>

#include "ckalloc.H"
#include "sim4parameters.H"
#include "sim4command.H"
#include "exon.H"
#include "mspManager.H"

#include "sim4defines.H"

//
//  A wrapper for the sim4 subroutines
//

enum { INIT = 0, PERM = 1, TEMP = 2};
enum { FWD = 0, BWD = 1, BOTH = 2 };
enum { OK = 0, FREE_START = 1, FREE_END = 2, FREE_BOTH_ENDS = 3};

/* data structures */

typedef struct coordinates {
  int pos1;
  int pos2;
}  coords;

/* used only in the alignment stage */
typedef struct edit_script {
  char op_type;   /* SUB, INS, or DEL */
  int num;        /* Number of operations */
  struct edit_script *next;
} edit_script;

typedef struct edit_script_list {
  int   offset1, offset2;
  int   len1, len2; 
  int    score;
  struct edit_script *script;
  struct edit_script_list *next_script;
} edit_script_list;

typedef  void  *Pointer;

typedef  struct ValNode {
  Pointer      data;
  struct ValNode *next;
} *ValNodePtr;

typedef struct spliced {
  int xs, xe, ys, ye, score;
  char type;
  struct spliced *next;
} splice_t;

typedef struct sim4_stats {
  int    internal;
  int    icoverage;
  int    numberOfMatches;
  int    numberOfNs;
  int    orientation;
  int    percentID;
#if ABORT_EXPENSIVE
  bool   tooManyMSPs;
#endif
} sim4_stats_t;



#define GEN_LOG4_ENTRIES 45
#define CDNA_LOG4_ENTRIES 25

#ifndef FOR_BRI_TABLE

#define HASH_SIZE 32767        /* 2**15 - 1 */

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

#endif






class Sim4 {
private:
  mspManager        _mspManager;

  sim4parameters   *globalParams;

  char       *seq1;
  char       *seq2;

  int         len1;
  int         len2;

  int         encoding[256];

  coords      last_GT;
  coords      last_CT;
  coords      last_AG;
  coords      last_AC;

  //sim4_args_t rs;

  int         X;
  int         wordSize;
  int         mspThreshold1;
  int         mspThreshold2;

  int         G_score;
  int         C_score;
  Exon       *exon_list;

  int mask;

#ifndef FOR_BRI_TABLE
  struct hash_table  *hashtable;
  struct hash_table   phashtable;
  struct hash_table   thashtable;
#endif


  //  For high-accuracy switch
  //
  bool _accurateSequences;

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

  int snake(int k, int x, int endx, int endy)
      {
        int y;

        if (x<0)
          return x;
        y = x+k;
        while (x<endx && y<endy && seq1[x]==seq2[y]) {
          ++x;
          ++y;
        }
        return x;
      };


  int rsnake(int k, int x, int startx, int starty, int m) {
        int y;

        if (x>m) return x;

#ifdef DEBUG
        if ((startx<0) || (starty<0))
          (void)printf("TROUBLE!!! startx:  %5d,  starty:  %5d\n",startx, starty);
        if ((x>M) || (x+k>N))
          (void)printf("TROUBLE!!! x:  %5d,  y:  %5d\n",x,x+k);
#endif
        
        y = x+k;
        while (x>startx && y>starty && seq1[x-1]==seq2[y-1]) {
          --x;
          --y;
        }
        return x;
      };


  char *checkExonsForOverlaps(Exon *theExons);
  char *appendExons(char *outputstring, Exon *theExon);
  void  maskExonsFromGenomic(Exon *theExons, char *f, char *r, int l);

  char   *IDISPLAY(char *outputstring,
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
  char   *appendAlignments(char *outputstring,
                           char *seq1,
                           char *seq2,
                           int len1,
                           int len2, 
                           edit_script_list **Aligns,
                           Exon *Exons, 
                           int match_ori);






  //  log.C
  int    getMSPthreshold(int defaultK, int len1, int len2);


  //  Functions from sim4b1.h
  //
  void link_to_data_list(Pointer,ValNodePtr *,ValNodePtr *);
  void ValNodeFreeData(ValNodePtr);
  int good_ratio(int);
  void flip_list(Exon **,Exon **);
  void free_list(Exon *);
  void free_align(edit_script_list *);
  void complement_exons(Exon **,int,int);
#if 0
  void print_exons(FILE *, Exon *);
#endif


#ifndef FOR_BRI_TABLE
  void bld_table(char *,int,int,int);
#else
  void bld_table(char *,int,int,int) {}
#endif

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

  struct edit_script_list *SIM4(char*,char*,int,int,int*,Exon**,
                                int *pA,
                                int *pT,
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
  //int    find_log_entry(long const *log4s, int n, int len, int offset);
  //Exon  *new_exon(int,int,int,int,int,int,int,Exon *);
  void   add_word(int,int);
  int     extend_hit(int,int,const char *const,const char * const,int,int,int,int);
  void   search(char *,char *,int,int,int,int); 
  void   exon_cores(char*,char*,int,int,int,int,int,int,int,int);
  Exon  *find_previous(Exon *,Exon *);
  bool   get_match_quality(Exon *,Exon *,sim4_stats_t *,int);
  void   check_consistency_intron_ori(Exon *,int,char *);


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

#if ABORT_EXPENSIVE
    _mspManager.setLimits(p->_mspLimitAbsolute, p->_mspLimitPercent);
#endif
    _mspManager.setParameters(p->_match, p->_mismatch, p->_percentError);


    seq1 = 0L;
    seq2 = 0L;

    wordSize      = DEFAULT_W;
    mspThreshold1 = p->_mspThresh1;  //K;
    mspThreshold2 = p->_mspThresh2;  //C;

    X        = DEFAULT_X;

    G_score = 0;
    C_score = 0;

    mask = 0;

#ifndef FOR_BRI_TABLE
    hashtable = 0L;
#endif

#if 0
    for (u32bit i=0; i<HASH_SIZE+1; i++) {
      phashtable.table[i] = 0L;
      thashtable.table[i] = 0L;
    }
#endif

#ifndef FOR_BRI_TABLE
    phashtable.nextPos   = 0L;
    phashtable.nodes     = 0L;
    phashtable.nodesused = 0;

    thashtable.nextPos   = 0L;
    thashtable.nodes     = 0L;
    thashtable.nodesused = 0;
#endif

    //  Initialize the encoding
    //
    for (u32bit i=256; i;)
      encoding[--i] = -1;

    encoding['A'] = encoding['a'] = 0;
    encoding['C'] = encoding['c'] = 1;
    encoding['G'] = encoding['g'] = 2;
    encoding['T'] = encoding['t'] = 3;

    _accurateSequences = false;

    last_AG.pos1 = last_AG.pos2 = 0;
    last_AC.pos1 = last_AC.pos2 = 0;
    last_GT.pos1 = last_GT.pos2 = 0;
    last_CT.pos1 = last_CT.pos2 = 0;
  };

  ~Sim4() {
#ifndef FOR_BRI_TABLE
    delete [] phashtable.nextPos;
    delete [] phashtable.nodes;

    delete [] thashtable.nextPos;
    delete [] thashtable.nodes;
#endif
  };

  char  *run(sim4command *cmd);
};


#endif  //  SIM4_H
