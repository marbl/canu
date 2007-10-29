#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include "bio.h"
#include "sim4.H"

//  Kaz Kylheku <kaz@ashi.footprints.net> library.
#include "kazlib/dict.h"
#include "kazlib/except.h"
#include "kazlib/hash.h"
#include "kazlib/list.h"
#include "kazlib/sfx.h"

//  Derived from pickBestPolish.c.  We report only the single best
//  match, when it is obvious that there is EXACTLY one best match.
//
//  Example: we have ten matches, but one is 3%id better than everyone
//  else -- that is an obviously unique match.  The rest are noise.
//
//  Example: ten matches, but they're all about the same quality -- within
//  a few percent id, and about the same length.  We pick no match, and
//  silently discard all.
//

//  Modified to:
//    a)  not print out unique matches
//    b)  print hangs
//    c)  print q20 bases inside mapped regions, outside, etc.
//
//  It needs two args -f seq.fasta -q qlt.fasta, both must have an
//  index -- build it for the seq.fasta, and COPY the index to
//  qlt.fastaidx.  Be sure to 'touch -r seq.fasta qlt.fasta' to get
//  the same timestamp on the files.
//
//
//  Further modified to behave like pickUniquePolish (print unique matches
//  to a specific file).
//
//  so:  pickUniquePolish-nhgri needs to read polishes on stdin
//    -f qry.fasta       -- query sequences for quality comparison
//    -q qlt.fasta       -- 
//    -scores X.scores   -- write stats to file X
//    -unique X.bz2      -- write uniquely mapped stuff to bzip2 file X.bz2
//    -filter X          -- filter out polishes less than X% of the longest
//    -output X.bz2      -- write filtered polishes to bzip2 file X.bz2
//  
//  It has two modes:
//    -f -q         -- just compute stats on the input.
//    all options   -- filter, and compute stats.
//
//  bzip2 -dc pass?/map-gen*-qlt$id.sim4db.bz2 |
//  $bin/fixPolishesIID -c $qry -g $gen |
//  $bin/filterPolishes -node -D |
//  $bin/sortPolishes -c -m 768 -t /scratch -v |
//  $bin/pickUniquePolish-nhgri > all-$id.scores
//    -o all-$id.sim4db.bz2
//    -F X
//    -f $qry
//    -q $qlt
//    -stats   all-$id.scores |
//    -uniq    all-$id.sim4db.bz2



u32bit  statOneMatch      = 0;
u32bit  statConsistent    = 0;
u32bit  statInconsistent  = 0;
u32bit  statUnique        = 0;
u32bit  statLost          = 0;

u32bit  consistentTie         = 0;
u32bit  consistentMatches     = 0;
u32bit  consistentIdentity    = 0;
u32bit  consistentTooShort    = 0;
u32bit  consistentNot         = 0;

u32bit  totLQ = 0;
u32bit  totMQ = 0;
u32bit  totRQ = 0;

seqFile      *SEQ = 0L;
seqFile      *QLT = 0L;

double       filter      = 0.0;
FILE        *oFile       = 0L;
int          oFileIsPipe = 0;
FILE        *sFile       = 0L;
FILE        *uFile       = 0L;
bool         doFiltering = false;


void
analyze(u32bit   iid,
        u32bit   clrl,
        u32bit   clrr,
        u32bit   len,
        bool     isForward,
        char     type) {

  QLT->find(iid);
  seqInCore  *Q = QLT->getSequenceInCore();;

  char  *q = Q->sequence();

  u32bit i = 0;

  u32bit  lq = 0;
  u32bit  mq = 0;
  u32bit  rq = 0;

  for ( ;i<clrl; i++)
    if (q[i] >= '0' + 20)
      lq++;

  for ( ;i<clrr; i++)
    if (q[i] >= '0' + 20)
      mq++;

  for ( ; i<len; i++)
    if (q[i] >= '0' + 20)
      rq++;

  delete Q;

  if (isForward) {
    totLQ += lq;
    totMQ += mq;
    totRQ += rq;
  } else {
    totLQ += rq;
    totMQ += mq;
    totRQ += lq;
  }

  fprintf(sFile, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%c\n",
          iid, clrl, clrr, len, lq, mq, rq, type);
}


void
analyze(sim4polish *p,
        char        type) {

  u32bit  clrl = p->exons[0].estFrom - 1;
  u32bit  clrr = p->exons[0].estTo   - 1;

  if (p->matchOrientation == SIM4_MATCH_COMPLEMENT) {
    clrl = p->estLen - (p->exons[0].estTo   - 1);
    clrr = p->estLen - (p->exons[0].estFrom - 1);
  }

  analyze(p->estID, clrl, clrr, p->estLen, p->matchOrientation != SIM4_MATCH_COMPLEMENT, type);
}



void
pickBestSlave(sim4polish **p, u32bit pNum) {
  u32bit        identitym = 0, nmatchesm = 0;  //  Best score for the mList
  u32bit        identityi = 0, nmatchesi = 0;  //  Best score the the iList
  u32bit        matchi = 0,    matchm = 0;

  //  Difficult choice here....
  //
  if (pNum == 1) {
    statOneMatch++;
    statUnique++;

    if (uFile)
      s4p_printPolish(uFile, p[0], S4P_PRINTPOLISH_FULL);

    if (oFile)
      s4p_printPolish(oFile, p[0], S4P_PRINTPOLISH_FULL);

    analyze(p[0], 'U');

    return;
  }

  //  Find the best percentIdentity and best numberOfMatches.  
  //
  //  identityi is the best percent identity of all the matches for this EST, and
  //  nmatchesi is the number of matches for the longest best identity match(es).
  //  matchi    is the match index
  //
  //  nmatchesm is the best numMatches of all the matches for this EST, and 
  //  identitym is the highest percent identity for the best numMatches match(es).
  //  matchm    is the match index

  for (u32bit i=0; i<pNum; i++) {
    if ((p[i]->percentIdentity > identityi) || 
        (p[i]->percentIdentity == identityi && p[i]->numMatches > nmatchesi)) {
      identityi = p[i]->percentIdentity;
      nmatchesi = p[i]->numMatches;
      matchi    = i;
    }
   
    if ((p[i]->numMatches > nmatchesm) ||
        (p[i]->numMatches == nmatchesm && p[i]->percentIdentity > identitym)) {
      nmatchesm = p[i]->numMatches;
      identitym = p[i]->percentIdentity;
      matchm    = i;
    }
  }


  bool  matchIsOK = false;

  //  If we are in agreement on what the best quality match is,
  //  see if the best match is obviously unique.
  //
  if ((identityi == identitym) ||
      (nmatchesi == nmatchesm)) {
    statConsistent++;

    //  It's clear what the quality values of the best match is, but we
    //  don't know if those values are shared by more than one match.
    //  Count the number of matches with exactly those scores.  If
    //  there is more than one, then we cannot pick out a single best.
    //
    u32bit numBest = 0;
    for (u32bit i=0; i<pNum; i++)
      if ((p[i]->percentIdentity == identityi) && (p[i]->numMatches == nmatchesi))
        numBest++;

    if (numBest > 1) {

      //  Dang, we mapped this guy more than once, exactly the same!
      //
      consistentTie++;

    } else {

      //  We claim to have a single best match.  See if any other
      //  matches are close to the quality of that one.

      u32bit  closeQuality = 0;

      for (u32bit i=0; i<pNum; i++)
        if (((p[i]->percentIdentity * 102) >= (identityi * 100)) ||
            ((p[i]->numMatches      * 102) >= (nmatchesi * 100)))
          closeQuality++;

      //  If only one match has close quality (the one we want to save!),
      //  save it.  Otherwise, label this query as multiple.

      u32bit  length = p[matchi]->exons[0].estFrom - p[matchi]->exons[0].estTo;

      if (closeQuality == 1) {
        matchIsOK = true;
        consistentMatches++;
      } else if ((length > 100) &&
                 (length / p[matchi]->estLen < 0.5)) {
        consistentTooShort++;
      } else {
        consistentNot++;
      }
    }

  } else {

    //  Otherwise, we disagree on what the best match is.
    //
    //  That is, the match with the highest identity is not the match
    //  with the highest number of matches -- a longer match exists, but
    //  at lower overall percent identity.

    statInconsistent++;

    //  Estimate the identity of the extended part, assuming the piece
    //  matched in common is matched at about the same identity.  Or
    //  just give up and say it's mapped to multiple places!

  }


  u32bit  best  = 0;
  u32bit  besti = 0;
  
  if (matchIsOK) {
    statUnique++;
    if (uFile)
      s4p_printPolish(uFile, p[matchi], S4P_PRINTPOLISH_FULL);

    assert(matchi == matchm);

    besti = matchi;
    analyze(p[besti], 'G');
  } else {
    statLost++;

    //  Just pick the longest match, analyze that.

    for (u32bit i=0; i<pNum; i++) {
      u32bit  len = p[i]->exons[0].estFrom - p[i]->exons[0].estTo;

      if ((len  > best) ||
          ((len == best) && (p[i]->numMatches > p[besti]->numMatches))) {
        best  = len;
        besti = i;
      }
    }

    analyze(p[besti], 'N');
  }


#if 0
  u32bit  nm = (u32bit)(p[besti]->numMatches * 0.75);
  u32bit  sv = 0;

  for (u32bit i=0; i<pNum; i++)
    if (p[i]->numMatches >= nm)
      sv++;
  
  fprintf(stderr, "Saved "u32bitFMT" matches more than nmatches "u32bitFMT" (from best of "u32bitFMT")\n", sv, nm, p[besti]->numMatches);
#endif


  //  besti is the best/longest match we have.  Decide on a threshold
  //  to throw out the obvious junk.
  //
  if ((oFile) && (doFiltering)) {
    u32bit  nm = (u32bit)(p[besti]->numMatches * filter);

    for (u32bit i=0; i<pNum; i++)
      if (p[i]->numMatches >= nm)
        s4p_printPolish(oFile, p[i], S4P_PRINTPOLISH_FULL);
  }

#if 0
  fprintf(stderr, "Uni:"u32bitFMTW(8)" Con:"u32bitFMTW(8)" (T:"u32bitFMTW(8)" M:"u32bitFMTW(8)" I:"u32bitFMTW(8)" N:"u32bitFMTW(8)") Inc:"u32bitFMTW(8)" -- Save:"u32bitFMTW(8)" Lost:"u32bitFMTW(8)"\r",
          statOneMatch,
          statConsistent, consistentTie, consistentMatches, consistentIdentity, consistentNot,
          statInconsistent,
          statUnique, statLost);
#endif
}







//  Just a wrapper around the real best picker, so that we can easily
//  destroy polishes when we're done.
//
void
pickBest(sim4polish **p, u32bit pNum) {

  pickBestSlave(p, pNum);

  for (u32bit i=0; i<pNum; i++)
    s4p_destroyPolish(p[i]);
}







dict_t  *IIDdict = 0L;
dict_t  *SEQdict = 0L;
dict_t  *GENdict = 0L;

void
fixIID(sim4polish *q, dict_t *estdict) {

  //  Fix the IID's
  dnode_t *cid = dict_lookup(estdict, q->estDefLine);
  dnode_t *gid = dict_lookup(GENdict, q->genDefLine);

  if ((cid == 0L) || (gid == 0L)) {
    const char *msg = "both deflines";
    if (cid)  msg = "genomic defline";
    if (gid)  msg = "est defline";

    s4p_printPolish(stdout, q, S4P_PRINTPOLISH_NOTVALUABLE);
    fprintf(stderr, "ERROR:  Couldn't find %s (%p %p) in the dictionary!\n", msg, cid, gid);
    exit(1);
  }

  q->estID = (u32bit)dnode_get(cid);
  q->genID = (u32bit)dnode_get(gid);
}











//
//  Stolen from sortPolishes
//
int          mergeFilesLen;
int          mergeFilesMax;
FILE       **mergeFiles;
char       **mergeNames;
sim4polish **mergePolishes;

sim4polish *
nextPolish(void) {
  int smallestPolish = 0;
  int nextPolish     = 1;

  //  If no merge files, read from stdin
  //
  if (mergeFilesLen == 0) {
    return(s4p_readPolish(stdin));
  }

  //  Find the smallest polish.
  //
  for (nextPolish = smallestPolish+1; nextPolish < mergeFilesLen; nextPolish++) {
    if (s4p_estIDcompare(mergePolishes+smallestPolish, mergePolishes+nextPolish) > 0)
      smallestPolish = nextPolish;
  }

  //  If the smallestPolish is 0L, we're all done.  Otherwise, dump
  //  the current smallest and fill it with a new polish.
  //
  if (mergePolishes[smallestPolish] == 0L) {
    return(0L);
  } else {
    sim4polish  *ret = mergePolishes[smallestPolish];
    mergePolishes[smallestPolish] = s4p_readPolish(mergeFiles[smallestPolish]);

    //  fix the iid's to be consistent in our partition, so we can have the input files
    //  sorted by est iid.
    if (mergePolishes[smallestPolish])
      fixIID(mergePolishes[smallestPolish], IIDdict);

    //  fix the iid's to be consistent globally
    fixIID(ret, SEQdict);

    return(ret);
  }
}







//
//  Stolen from fixPolishesIID
//
void
addToDict(dict_t *d, char *n) {
  dnode_t  *node = 0L;
  char     *dcpy = 0L;

  if (n == 0L)
    return;

  seqFile  *F = openSeqFile(n);
  while (!F->eof()) {
    seqInCore *S = F->getSequenceInCore();

    node = (dnode_t *)palloc(sizeof(dnode_t));
    dcpy = (char    *)palloc(sizeof(char) * S->headerLength() + 1);

    strcpy(dcpy, S->header());

    dnode_init(node, (void *)S->getIID());
    dict_insert(d, node, dcpy);

    delete S;
  }
  delete F;
}

int
headerCompare(const void *a, const void *b) {
  char  *A = *((char **)a);
  char  *B = *((char **)b);

  //fprintf(stderr, "%s -- %s\n", A, B);
  return(strcmp(A, B));
}








int
main(int argc, char **argv) {
  u32bit       pNum   = 0;
  u32bit       pAlloc = 8388608;
  sim4polish **p      = 0L;
  sim4polish  *q      = 0L;
  u32bit       estID  = ~0;

  bool        *found  = 0L;

  //  From fixPolishesIID.c
  IIDdict = 0L;
  SEQdict = 0L;
  GENdict = 0L;

  //  Incorporated from sortPolishes
  mergeFilesLen   = 0;
  mergeFilesMax   = sysconf(_SC_OPEN_MAX);
  mergeFiles      = (FILE **)malloc(sizeof(FILE*) * mergeFilesMax);
  mergeNames      = (char **)malloc(sizeof(char*) * mergeFilesMax);
  mergePolishes   = (sim4polish **)malloc(sizeof(char*) * mergeFilesMax);

  //  Default to printing stats on stdout.
  sFile = stdout;

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-n") == 0) {
      pAlloc = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-fpart") == 0) {
      arg++;
      fprintf(stderr, "reading query deflines from '%s'\n", argv[arg]);
      IIDdict = dict_create(DICTCOUNT_T_MAX, headerCompare);
      addToDict(IIDdict, argv[arg]);
    } else if (strcmp(argv[arg], "-g") == 0) {
      ++arg;
      fprintf(stderr, "reading genomic deflines from '%s'\n", argv[arg]);
      GENdict = dict_create(DICTCOUNT_T_MAX, headerCompare);
      addToDict(GENdict, argv[arg]);
    } else if (strcmp(argv[arg], "-F") == 0) {
      ++arg;
      fprintf(stderr, "reading query deflines from '%s'\n", argv[arg]);
      SEQdict = dict_create(DICTCOUNT_T_MAX, headerCompare);
      addToDict(SEQdict, argv[arg]);
    } else if (strcmp(argv[arg], "-f") == 0) {
      ++arg;
      SEQ = openSeqFile(argv[arg]);
      SEQ->openIndex();
    } else if (strcmp(argv[arg], "-q") == 0) {
      ++arg;
      QLT = openSeqFile(argv[arg]);
      QLT->openIndex();

    } else if (strcmp(argv[arg], "-filter") == 0) {
      filter = atof(argv[++arg]);
      doFiltering = true;
    } else if (strcmp(argv[arg], "-output") == 0) {
      char  cmd[1024] = {0};
      errno = 0;
      ++arg;
      if (strcmp(argv[arg] + strlen(argv[arg]) - 4, ".bz2") == 0) {
        sprintf(cmd, "bzip2 -1c > %s", argv[arg]);
        oFile = popen(cmd, "w");
        oFileIsPipe = 1;
      } else if (strcmp(argv[arg] + strlen(argv[arg]) - 3, ".gz") == 0) {
        sprintf(cmd, "gzip -1c > %s", argv[arg]);
        oFile = popen(cmd, "w");
        oFileIsPipe = 1;
      } else {
        fprintf(stderr, "Got %s, not .bz2 not .gz!\n", argv[arg]);
        exit(1);
      }
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", cmd, strerror(errno));
      doFiltering = true;
    } else if (strcmp(argv[arg], "-scores") == 0) {
      errno = 0;
      sFile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", argv[arg-1], strerror(errno));
      doFiltering = true;
    } else if (strcmp(argv[arg], "-unique") == 0) {
      char  cmd[1024] = {0};
      errno = 0;
      arg++;
      if (strcmp(argv[arg] + strlen(argv[arg]) - 4, ".bz2") == 0)
        sprintf(cmd, "bzip2 -1c > %s", argv[arg]);
      else if (strcmp(argv[arg] + strlen(argv[arg]) - 3, ".gz") == 0)
        sprintf(cmd, "gzip -1c > %s", argv[arg]);
      else
        sprintf(cmd, "cat > %s", argv[arg]);
      uFile = popen(cmd, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", cmd, strerror(errno));
      doFiltering = true;

    } else if (strncmp(argv[arg], "-M", 2) == 0) {
      arg++;
      while ((arg < argc) && (fileExists(argv[arg]))) {
        if (mergeFilesLen >= mergeFilesMax) {
          fprintf(stderr, "%s: ERROR!  Too many input files!  Should be less than %d\n", argv[0], mergeFilesMax);
          exit(1);
        }
        mergeNames[mergeFilesLen]   = argv[arg];
        mergeFiles[mergeFilesLen++] = openFile(argv[arg], "r");
        arg++;
      }
      arg--;

    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }


  if (doFiltering) {
    if (uFile == 0L)
      fprintf(stderr, "ERROR:  -unique is required\n"), exit(1);
    if (sFile == 0L)
      fprintf(stderr, "ERROR:  -scores is required\n"), exit(1);
    if ((filter < 0.0) || (filter > 1.0))
      fprintf(stderr, "ERROR:  -filter value of %f invalid.  0 <= F <= 100.\n", filter), exit(1);
  }


  if ((IIDdict == 0L) || (SEQdict == 0L) || (GENdict == 0L)) {
    fprintf(stderr, "WARNING!  No sequence dictionaries, NOT FIXING IIDs!  (supply -fpart, -f and -g)\n");
  }


  if ((SEQ == 0L) || (QLT == 0L)) {
    fprintf(stderr, "I need -f and -q\n");
    exit(1);
  }

  //  We no longer require that input polishes be sorted increasingly;
  //  now they only must be grouped.  This remembers if we've seen a
  //  match or not.  At the end, we'll analyze() those we haven't done
  //  already.
  //
  found = new bool [ SEQ->getNumberOfSequences() ];
  for (u32bit i=0; i<SEQ->getNumberOfSequences(); i++)
    found[i] = false;


  //  Initialize the merge -- if no merge files, nothing done!
  //
  for (int i=0; i<mergeFilesLen; i++) {
    mergePolishes[i] = s4p_readPolish(mergeFiles[i]);
    fixIID(mergePolishes[i], IIDdict);
  }


  //  Read polishes, picking the best when we see a change in the
  //  estID.

  p = (sim4polish **)malloc(sizeof(sim4polish *) * pAlloc);

  while ((q = nextPolish()) != 0L) {

    if ((q->estID != estID) && (pNum > 0)) {
      //fprintf(stderr, "PickBest for estID "u32bitFMT"\n", estID);

      found[estID] = true;
      pickBest(p, pNum);
      pNum  = 0;
    }

    //  Reallocate pointers?
    //
    if (pNum >= pAlloc) {
      p = (sim4polish **)realloc(p, sizeof(sim4polish *) * (pAlloc *= 2));
      if (p == 0L) {
        fprintf(stderr, "Out of memory: Couldn't allocate space for polish pointers.\n");
        exit(1);
      }
    }

    p[pNum++] = q;
    estID     = q->estID;
  }

  if (pNum > 0) {
    found[estID] = true;
    pickBest(p, pNum);
  }

  //  Attempt cleanup
  //
  for (int i=0; i<mergeFilesLen; i++)
    closeFile(mergeFiles[i], mergeNames[i]);

  for (estID=0; estID < SEQ->getNumberOfSequences(); estID++)
    if (found[estID] == false)
      analyze(estID, 0, SEQ->sequenceLength(estID), SEQ->sequenceLength(estID), true, 'M');

  if (oFile)  pclose(oFile);
  if (uFile)  pclose(uFile);
  if (sFile)  fclose(sFile);

  fprintf(stderr, "Uni:"u32bitFMTW(8)" Con:"u32bitFMTW(8)" (T:"u32bitFMTW(8)" M:"u32bitFMTW(8)" I:"u32bitFMTW(8)" S:"u32bitFMTW(8)" N:"u32bitFMTW(8)") Inc:"u32bitFMTW(8)" -- Save:"u32bitFMTW(8)" Lost:"u32bitFMTW(8)"\n",
          statOneMatch,
          statConsistent, consistentTie, consistentMatches, consistentIdentity, consistentTooShort, consistentNot,
          statInconsistent,
          statUnique, statLost);
  fprintf(stderr, "total:  LQ:"u32bitFMT" MQ:"u32bitFMT" RQ:"u32bitFMT"\n",
          totLQ, totMQ, totRQ);

  return(0);
}

