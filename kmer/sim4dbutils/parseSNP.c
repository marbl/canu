#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#include "bio.h"
#include "sim4.H"

//  Writes things with mappings that don't contain the snp itself to a
//  failure file.  Otherwise, if the mapping is above the threshold, a
//  line describing the snp is output.


char *usage =
"usage: %s [options]\n"
"             -i min-identity     filter matches on percent identity\n"
"             -c min-coverage     filter matches on percent coverage\n"
"             -F failed           save matches that do not contain the\n"
"                                 to the file 'failed'\n"
"             -O output           save the parsed SNPs to the file\n"
"                                 'output'\n"
"             -D prefix           report debugging stuff into files\n"
"                                 prefixed with 'prefix'\n"
"             -d delimiter        Use the single character delimiter as\n"
"                                 the end of the defline ID field.  The\n"
"                                 default is to split on any whitespace.\n"
"\n"
"             -s sizeTag          Use this tag as the size of the snp.\n"
"                                 '/size=' is tried by default.\n"
"\n"
"             -p posTag           Use this tag as the position of the snp.\n"
"                                 'allelePos=' and '/pos=' are tried by\n"
"                                 default, and if posTag is not found.\n"
"\n"
"                                 TAGS: The number immediately after the first\n"
"                                 occurance of the tag will be used.\n"
"\n"
"             -o offset           An additive offset to the SNP position.\n"
"                                 The default is 0.\n"
"\n"
"             -format n           1 - use the original (default) format\n"
"                                 2 - use an extended format, includes the\n"
"                                     position in the alignment string\n"
"\n"
"             -h                  Show this help.\n"
"\n"
"\n" 
"             only -O is required.  Input is read from stdin.\n"
"\n"
"             NOTE!  Sizes and sizeTag is NOT IMPLEMENTED!\n"
"                    All SNPs are of size == 1\n"
"\n"
"             If you parse base-based SNPs, the result is returned base-based.\n"
"             You should use an ofset of 0.\n"
"\n"
"             If you parse space-based SNPs, the result is returned base-based.\n"
"             You should use an offset of 1.\n"
"\n";

FILE *multiMultiFile   = 0L;  //  multiple hits, at least one is multiple exon
FILE *multiSingleFile  = 0L;  //  multiple hits, all are single exon
FILE *singleMultiFile  = 0L;  //  single hit, and it has more than one exon
FILE *singleSingleFile = 0L;  //  single hit, single exon

int   smpass = 0;
int   sspass = 0;
int   mmpass = 0;
int   mspass = 0;

int   smfail = 0;
int   ssfail = 0;
int   mmfail = 0;
int   msfail = 0;

int   failedsnps    = 0;
int   failedmatches = 0;

FILE *validSNPMapFile   = 0L;
FILE *failedSNPMapFile  = 0L;

char  fieldDelimiter = 0;
char *sizeTag        = "/size=";
char *posTag         = "/pos=";
int   positionOffset = 0;

int   outputFormat   = 1;

char *
findSNPid(char *defline) {
  char *ret = 0L;
  int   sta = 0;
  int   len = 0;
  int   i = 0;

  if (fieldDelimiter == 0) {
    for (len=1; defline[len] && !isspace(defline[len]); len++)
      ;
  } else {
    for (len=1; defline[len] && defline[len] != fieldDelimiter; len++)
      ;
  }

#if 0
  //  This was used for a set of SNPs with a non-standard defline
  //  structure.  It returns the field between the first '|' and the
  //  next '_'.
  //
  for (len=1; defline[len] && defline[len] != '_'; len++)
    ;
  for (sta=len-1; sta > 0 && defline[sta] != '|'; sta--)
    ;
#endif

  errno = 0;
  ret = (char *)malloc(sizeof(char) * (len+1));
  if (errno) {
    fprintf(stderr, "malloc() problem: %s\n", strerror(errno));
    exit(1);
  }

  for (i=sta; i<len-1; i++)
    ret[i-sta] = defline[i+1];

  ret[len-sta-1] = 0;

  return(ret);
}

char *
findGENid(char *defline) {
  char *ret = 0L;
  int   len = 0;
  int   i = 0;

  for (len=1; defline[len] && !isspace(defline[len]); len++)
    ;

  errno = 0;
  ret = (char *)malloc(sizeof(char) * (len+1));
  if (errno) {
    fprintf(stderr, "malloc() problem: %s\n", strerror(errno));
    exit(1);
  }

  for (i=0; i<len-1; i++)
    ret[i] = defline[i+1];

  ret[len-1] = 0;

  return(ret);
}



int
findPosition(char *defline) {
  char *p = 0L;

  p = strstr(defline, posTag);

  //  Look for standard posTags if we didn't find the one the user wanted.

  if (p == 0L)
    p = strstr(defline, "allelePos=");
  if (p == 0L)
    p = strstr(defline, "/pos=");
  if (p == 0L)
    fprintf(stderr, "posTag '%s' (also looked for 'allelePos=' and '/pos=') not found in defline '%s'!\n", posTag, defline), exit(1);

  while (*p && !isdigit(*p))
    p++;

  if (*p == 0) {
    fprintf(stderr, "Found posTag '%s' in defline '%s', but didn't find any numbers!\n", posTag, defline);
    exit(1);
  }

  return(atoi(p) + positionOffset);
}




//  Returns 1 if SNP was valid and printed,
//  0 otherwise.
//
int
printSNP(FILE *F, sim4polish *p) {
  u32bit   pos           = findPosition(p->estDefLine);
  u32bit   exonWithSNP   = ~u32bitZERO;
  u32bit   i             = 0;
  u32bit   seqOffset     = 0;

  //  If the match is complement, then the alignment is printed using
  //  the reverse complemented SNP sequence, and so we need to find
  //  the offset at the end of the sequence (not always the same as
  //  the offset at the start of the sequence).
  //
  //  XXX: Previous version had this as "p->estLen - pos + siz", which
  //  seems wrong.  This version does what appears to be reverse
  //  complement - size.  I don't understand if this is a "size" or
  //  just a "1" thing.
  //
  seqOffset = pos;
  if (p->matchOrientation == SIM4_MATCH_COMPLEMENT)
    seqOffset = p->estLen - pos - 1;

  //  Find the exon with the SNP
  //
  for (i=0; i<p->numExons; i++)
    if (((p->exons[i].estFrom-1) <= seqOffset) && (seqOffset <= (p->exons[i].estTo-1)))
      exonWithSNP = i;

  if (exonWithSNP == ~u32bitZERO)
    return(0);

  //  If we are printing to a file, continue to find the location, otherwise,
  //  just return.
  //
  if (F) {
    char *SNPid = findSNPid(p->estDefLine);
    char *GENid = findGENid(p->genDefLine);

    char  SNPbase = 0;
    char  GENbase = 0;

    //  Now, we examine the alignment strings to decide exactly
    //  where the SNP is located in the genomic.
    //
    //  bpToExaine - the number of bases we need to skip in the
    //  alignment (counted in the snp), +1 because we are currently at
    //  the bp before the alignment (so we need to skip one more space).
    //
    //  XXX:  these used to be int!
    //
    u32bit  bpToExamine = seqOffset - (p->exons[exonWithSNP].estFrom - 1) + 1;
    u32bit  examinePos  = 0;
    u32bit  genPosition = p->genLo + p->exons[exonWithSNP].genFrom - 1;

    //  Recent runs of dbSNP showed that we are off by one (too many if forward, too few if complement).  This is a hack to fix it.
    //
    if (p->matchOrientation == SIM4_MATCH_COMPLEMENT)
      bpToExamine++;
    else
      bpToExamine--;

    while (bpToExamine > 0) {

      //  If the SNP alignment eats up a base pair, decrement
      //  the number of bp left to examine.
      //
      if (p->exons[exonWithSNP].estAlignment[examinePos] != '-')
        bpToExamine--;

      //  If the the genomic alignment is not a gap, increment the
      //  position.
      //
      if (p->exons[exonWithSNP].genAlignment[examinePos] != '-')
        genPosition++;

      examinePos++;
    }

    //  Adjust the quality values, treating the SNP as a match always.
    //
    SNPbase = p->exons[exonWithSNP].estAlignment[examinePos-1];
    GENbase = p->exons[exonWithSNP].genAlignment[examinePos-1];

    p->exons[exonWithSNP].estAlignment[examinePos-1] = 'A';
    p->exons[exonWithSNP].genAlignment[examinePos-1] = 'A';

    s4p_updateAlignmentScores(p);

    p->exons[exonWithSNP].estAlignment[examinePos-1] = SNPbase;
    p->exons[exonWithSNP].genAlignment[examinePos-1] = GENbase;


    if (outputFormat == 1) {
      fprintf(F, "%s %s "u32bitFMT" %c/%c %s global["u32bitFMT" "u32bitFMT"] exon["u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"]\n",
              SNPid,
              GENid,
              genPosition,
              SNPbase,
              GENbase,
              (p->matchOrientation == SIM4_MATCH_FORWARD) ? "forward" : "complement",
              p->percentIdentity,
              p->querySeqIdentity,
              p->numExons,
              exonWithSNP,
              p->exons[exonWithSNP].percentIdentity,
              (u32bit)floor(100.0 * (double)p->exons[exonWithSNP].numMatches / (double)p->estLen));
    } else if (outputFormat == 2) {

      //  The format is all on one line, data fields separated by tab.
      //  No spaces -- "sa=C" instead of "sa = C"
      //
      //    SNPid
      //    GENid
      //    genomic position of SNP
      //    sa=c        -- snp allele
      //    ga=c        -- genome allele
      //    mo={f|r}    -- mapping orientation
      //    pi=n        -- percent identity
      //    pc=n        -- percent coverage
      //    nb=n        -- number of alignment blocks
      //    bl=n        -- alignment block with the snp
      //    bp=n        -- position of the snp in the alignment block
      //    bi=n        -- percent identity of the block
      //    bc=n        -- percent coverage of the block
      //
      //  The first three items are mandatory, are always in that
      //  order, and are always the first three.  The others are
      //  optional, and can occur in any order.  There might be more
      //  present than listed here.
      //
      //  The order and content should be consistent for any given
      //  version of the software.
      //
      fprintf(F, "%s %s "u32bitFMT" sa=%c ga=%c mo=%c pi="u32bitFMT" pc="u32bitFMT" nb="u32bitFMT" bl="u32bitFMT" bp="u32bitFMT" bi="u32bitFMT" bc="u32bitFMT"\n",
              "a", //SNPid,
              "b", //GENid,
              genPosition,
              p->exons[exonWithSNP].estAlignment[examinePos-1],                                       // sa
              p->exons[exonWithSNP].genAlignment[examinePos-1],                                       // ga
              (p->matchOrientation == SIM4_MATCH_FORWARD) ? 'f' : 'r',                                // mo
              p->percentIdentity,                                                                     // pi
              p->querySeqIdentity,                                                                    // pc
              p->numExons,                                                                            // nb
              exonWithSNP,                                                                            // bl
              examinePos,                                                                             // bp
              p->exons[exonWithSNP].percentIdentity,                                                  // bi
              (u32bit)floor(100.0 * (double)p->exons[exonWithSNP].numMatches / (double)p->estLen));   // bc
    } else {
    }

    free(SNPid);
    free(GENid);
  }

  return(1);
}



//  Just a wrapper around the real best picker, so that we can easily
//  destroy polishes when we're done.
//
void
parseSNP(sim4polish **p, int pNum) {
  int   numMulti  = 0;
  int   i;

  //  Count the number of matches that have more than one exon
  //
  for (i=0; i<pNum; i++)
    if (p[i]->numExons > 1)
      numMulti++;

  if (pNum == 1) {

    //
    //  Exactly one match for this SNP
    //

    if (numMulti == 0) {

      //  Match has one exon

      if (singleSingleFile)
        s4p_printPolish(singleSingleFile, p[0], S4P_PRINTPOLISH_FULL);

      if (printSNP(validSNPMapFile, p[0])) {
        sspass++;
      } else {
        ssfail++;
        if (failedSNPMapFile)
          s4p_printPolish(failedSNPMapFile, p[0], S4P_PRINTPOLISH_FULL);
      }
    } else {

      //  Match has more than one exon

      if (singleMultiFile)
        s4p_printPolish(singleMultiFile, p[0], S4P_PRINTPOLISH_FULL);

      if (printSNP(validSNPMapFile, p[0])) {
        smpass++;
      } else {
        smfail++;
        if (failedSNPMapFile)
          s4p_printPolish(failedSNPMapFile, p[0], S4P_PRINTPOLISH_FULL);
      }
    }
  } else {

    //
    //  More than one match for this SNP
    //

    if (numMulti == 0) {
      int pass=0, fail=0;

      //  All the matches are single exon

      if (multiSingleFile)
        for (i=0; i<pNum; i++)
          s4p_printPolish(multiSingleFile, p[i], S4P_PRINTPOLISH_FULL);

      for (i=0; i<pNum; i++)
        if (printSNP(validSNPMapFile, p[i])) {
          pass++;
        } else {
          fail++;
          if (failedSNPMapFile)
            s4p_printPolish(failedSNPMapFile, p[i], S4P_PRINTPOLISH_FULL);
        }

      if (pass==1)       sspass++;
      if (pass > 1)      mspass++;
      if (!pass && fail) msfail++;
    } else {
      int pass=0, fail=0;

      //  At least one match has more than one exon -- the correct one
      //  might be a single exon, but we don't know which is which.

      if (multiMultiFile)
        for (i=0; i<pNum; i++)
          s4p_printPolish(multiMultiFile, p[i], S4P_PRINTPOLISH_FULL);

      for (i=0; i<pNum; i++)
        if (printSNP(validSNPMapFile, p[i])) {
          pass++;
        } else {
          fail++;
          if (failedSNPMapFile)
            s4p_printPolish(failedSNPMapFile, p[i], S4P_PRINTPOLISH_FULL);
        }

      if (pass==1)       smpass++;
      if (pass > 1)      mmpass++;
      if (!pass && fail) mmfail++;
    }
  }

  for (i=0; i<pNum; i++)
    s4p_destroyPolish(p[i]);
}


int
main(int argc, char **argv) {
  int          arg    = 1;
  int          pNum   = 0;
  int          pAlloc = 8388608;
  sim4polish **p      = 0L;
  sim4polish  *q      = 0L;
  u32bit       estID  = 0;

  int          percentID = 0;
  int          percentCO = 0;

  validSNPMapFile   = 0L;
  failedSNPMapFile  = 0L;

  while (arg < argc) {
    if        (strncmp(argv[arg], "-i", 2) == 0) {
      percentID = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-c", 2) == 0) {
      percentCO = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-F", 2) == 0) {
      errno = 0;
      failedSNPMapFile = fopen(argv[++arg], "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
    } else if (strncmp(argv[arg], "-O", 2) == 0) {
      errno = 0;
      validSNPMapFile = fopen(argv[++arg], "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
    } else if (strncmp(argv[arg], "-D", 2) == 0) {
      char   name[1025];

      errno = 0;

      arg++;

      sprintf(name, "%s-multi-multi", argv[arg]);
      multiMultiFile = fopen(name, "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", name, strerror(errno));
        exit(1);
      }

      sprintf(name, "%s-multi-single", argv[arg]);
      multiSingleFile = fopen(name, "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", name, strerror(errno));
        exit(1);
      }

      sprintf(name, "%s-single-multi", argv[arg]);
      singleMultiFile = fopen(name, "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", name, strerror(errno));
        exit(1);
      }

      sprintf(name, "%s-single-single", argv[arg]);
      singleSingleFile = fopen(name, "w");
      if (errno) {
        fprintf(stderr, "Couldn't open '%s' for writing.\n%s\n", name, strerror(errno));
        exit(1);
      }
    } else if (strncmp(argv[arg], "-d", 2) == 0) {
      arg++;
      fieldDelimiter = argv[arg][0];
    } else if (strncmp(argv[arg], "-p", 2) == 0) {
      arg++;
      posTag = argv[arg];
    } else if (strncmp(argv[arg], "-s", 2) == 0) {
      arg++;
      sizeTag = argv[arg];
    } else if (strncmp(argv[arg], "-o", 2) == 0) {
      arg++;
      positionOffset = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-format", 2) == 0) {
      arg++;
      outputFormat = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-h", 2) == 0) {
      fputs(usage, stderr);
      exit(1);
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }


  if ((outputFormat != 1) && (outputFormat != 2)) {
    fprintf(stderr, "Invalid output format.  Must be 1 or 2.\n");
    exit(1);
  }


  //  Read polishes, parsing when we see a change in the estID.
  //  Really, we could parse one by one, but it's nice to know if the
  //  thing mapped more than once.
  //
  //  We could also extend this to discard matches that look
  //  suspicious -- or maybe pick the single best match for each.

  p = (sim4polish **)malloc(sizeof(sim4polish *) * pAlloc);

  while ((q = s4p_readPolish(stdin)) != 0L) {
    if (q->estID < estID) {
      fprintf(stderr, "ERROR:  Polishes not sorted by SNP idx!  this="u32bitFMT", looking for "u32bitFMT"\n",
              q->estID, estID);
      exit(1);
    }

    if ((q->estID != estID) && (pNum > 0)) {
      parseSNP(p, pNum);
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

    estID     = q->estID;

    if ((q->percentIdentity >= percentID) &&
        (q->querySeqIdentity >= percentCO)) {
      p[pNum++] = q;
    } else {
      s4p_destroyPolish(q);
    }
  }

  if (pNum > 0)
    parseSNP(p, pNum);

  fprintf(stdout, "SNPs with:\n");
  fprintf(stdout, "  single hit, single exon:        %6d\n", sspass);
  fprintf(stdout, "  single hit, multiple exons:     %6d\n", smpass);
  fprintf(stdout, "  multiple hits, single exon:     %6d\n", mspass);
  fprintf(stdout, "  multiple hits, multiple exons:  %6d\n", mmpass);
  fprintf(stdout, "SNPs that failed:\n");
  fprintf(stdout, "  single hit, single exon:        %6d\n", ssfail);
  fprintf(stdout, "  single hit, multiple exons:     %6d\n", smfail);
  fprintf(stdout, "  multiple hits, single exon:     %6d\n", msfail);
  fprintf(stdout, "  multiple hits, multiple exons:  %6d\n", mmfail);

  fclose(validSNPMapFile);
  fclose(failedSNPMapFile);

  fclose(multiMultiFile);
  fclose(multiSingleFile);
  fclose(singleMultiFile);
  fclose(singleSingleFile);

  return(0);
}

