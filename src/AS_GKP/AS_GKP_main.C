
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

#include "AS_global.H"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#include <map>

#include "AS_PER_genericStore.H"
#include "AS_PER_gkpStore.H"
#include "AS_UTL_fileIO.H"
#include "AS_MSG_pmesg.H"
#include "AS_GKP_include.H"

char            *progName       = NULL;
long		srandSeed	= 0;

gkStore         *gkpStore       = NULL;
gkFragment      *gkFrag1        = NULL;
gkFragment      *gkFrag2        = NULL;
FILE            *errorFP        = NULL;

char             fastqUIDmapName[FILENAME_MAX];
FILE            *fastqUIDmap = NULL;

static
void
usage(char *filename, int longhelp) {
  fprintf(stdout, "usage1: %s -o gkpStore [append/create options] <input.frg> <input.frg> ...\n", filename);
  fprintf(stdout, "usage2: %s -P partitionfile gkpStore\n", filename);
  fprintf(stdout, "usage3: %s [id-selection] [options] [format] gkpStore\n", filename);
  fprintf(stdout, "\n");
  fprintf(stdout, "----------------------------------------------------------------------\n");
  fprintf(stdout, "The first usage will append to or create a GateKeeper store:\n");
  fprintf(stdout, "  -a                     append to existing store\n");
  fprintf(stdout, "  -o <gkpStore>          append to or create gkpStore\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  -T                     do not check minimum length (for OBT)\n");
  fprintf(stdout, "  -F                     fix invalid insert size estimates\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  -v <vector-info>       load vector clear ranges into each read.\n");
  fprintf(stdout, "                         MUST be done on an existing, complete store.\n");
  fprintf(stdout, "                         example: -a -v vectorfile -o that.gkpStore\n");
  fprintf(stdout, "                         format: 'UID vec-clr-begin vec-clr-end'\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "----------------------------------------------------------------------\n");
  fprintf(stdout, "The second usage will partition an existing store, allowing\n");
  fprintf(stdout, "the entire store partition to be loaded into memory.\n");
  fprintf(stdout, "  -P <partitionfile>     a list of (partition fragiid)\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "----------------------------------------------------------------------\n");
  fprintf(stdout, "The third usage will dump the contents of a GateKeeper store.\n");
  fprintf(stdout, "There are THREE components to a dump, what to dump, options, and format.\n");
  fprintf(stdout, "The first two are optional, the last is mandatory.  Examples:\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  Dump metainfo for the first 100 fragments\n");
  fprintf(stdout, "    gatekeeper -b 1 -e 100 -tabular -dumpfragments my.gkpStore > first100.tsv\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  Dump a random 25%% of the reads in the first library\n");
  fprintf(stdout, "    gatekeeper -randomsubset 1 0.25 -dumpfrg my.gkpStore > random25.frg\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  Dump fasta sequence for the UIDs in 'uidFile'\n");
  fprintf(stdout, "    gatekeeper -uid uidFile -dumpfastaseq file -dumpfrg my.gkpStore\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  -----------------------------------\n");
  fprintf(stdout, "  [selection of what objects to dump]\n");
  fprintf(stdout, "  -----------------------------------\n");
  fprintf(stdout, "  -b <begin-iid>          dump starting at this library or read\n");
  fprintf(stdout, "  -e <ending-iid>         dump stopping after this iid\n");
  fprintf(stdout, "  -uid <uid-file>         dump only objects listed in 'uid-file'\n");
  fprintf(stdout, "  -iid <iid-file>         dump only objects listed in 'iid-file'\n");
  fprintf(stdout, "  -randommated  <lib> <n> pick n mates (2n frags) at random from library lib\n");
  fprintf(stdout, "  -randomsubset <lib> <f> dump a random fraction f (0.0-1.0) of library lib\n");
  fprintf(stdout, "  -randomlength <lib> <l> dump a random fraction of library lib, fraction picked\n");
  fprintf(stdout, "                          so that the untrimmed length is close to l bp\n");
  fprintf(stdout, "  -longestovermin <lib> <n> pick all reads longer than n bp from library lib\n");
  fprintf(stdout, "  -longestlength  <lib> <n> pick longest reads from library lib to add up to n total bp\n"); 
  fprintf(stdout, "  -deterministic	     use a constant seed for random subset dumps\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  ---------\n");
  fprintf(stdout, "  [options]\n");
  fprintf(stdout, "  ---------\n");
  fprintf(stdout, "  -tabular               dump info, libraries or fragments in a tabular\n");
  fprintf(stdout, "                         format (for -dumpinfo, -dumplibraries,\n");
  fprintf(stdout, "                         and -dumpfragments, ignores -withsequence and -clear)\n");
  fprintf(stdout, "  -isfeatureset <libID> <X> sets exit value to 0 if feature X is set in library libID, 1 otherwise.\n");
  fprintf(stdout, "                        If libID == 0, check all libraries.\n");
  fprintf(stdout, "  -nouid                dump info without including the read UID (for -dumpinfo, -dumplibraries, -dumpfragments)\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  ----------------\n");
  fprintf(stdout, "  [format of dump]\n");
  fprintf(stdout, "  ----------------\n");
  fprintf(stdout, "  -dumpinfo                  print information on the store\n");
  fprintf(stdout, "    -lastfragiid             just print the last IID in the store\n");
  fprintf(stdout, "  -dumplibraries             dump all library records\n");
  fprintf(stdout, "  -dumpfragments             dump fragment info, no sequence\n");
  fprintf(stdout, "    -withsequence              ...and include sequence\n");
  fprintf(stdout, "    -clear <clr>               ...in clear range <clr>, default=LATEST\n");
  fprintf(stdout, "    -invert                    ...invert \n");
  fprintf(stdout, "  -dumpfasta <prefix>        dump fragment sequence and quality into <p.fasta> and <p.fasta.qual>\n");
  fprintf(stdout, "    -allreads                  ...all reads, regardless of deletion status (deleted are lowercase)\n");
  fprintf(stdout, "    -allbases                  ...all bases (lowercase for non-clear)\n");
  fprintf(stdout, "    -decoded                   ...quality as integers ('20 21 19')\n");
  fprintf(stdout, "    -clear <clr>               ...in clear range <clr>, default=LATEST\n");
  fprintf(stdout, "  -dumpfrg                   extract LIB, FRG and LKG messages\n");
  fprintf(stdout, "    -allreads                  ...all reads, regardless of deletion status\n");
  fprintf(stdout, "    -donotfixmates             ...only extract the fragments given, do not add in\n");
  fprintf(stdout, "                                  missing mated reads\n");
  fprintf(stdout, "    -clear <clr>               ...use clear range <clr>, default=LATEST\n");
  fprintf(stdout, "    -legacyformat              ...extract using frg format version 1 (legacy format, for compatibility)\n");
  fprintf(stdout, "  -dumpnewbler <prefix>      extract LIB, FRG and LKG messages, write in a\n");
  fprintf(stdout, "                             format appropriate for Newbler.  This will create\n");
  fprintf(stdout, "                             files 'prefix.fna' and 'prefix.fna.qual'.  Options\n");
  fprintf(stdout, "                             -donotfixmates and -clear also apply.\n");
  fprintf(stdout, "  -dumpfastq <prefix>        extract LIB, FRG and LKG messages, write in FastQ format.  Currently\n");
  fprintf(stdout, "                             this works only on a store with one library as all the mated reads\n");
  fprintf(stdout, "                             are dumped into a single file. This will create files 'prefix.paired.fastq',\n");
  fprintf(stdout, "                             'prefix.1.fastq', 'prefix.2.fastq' and 'prefix.unmated.fastq' for unmated\n");
  fprintf(stdout, "                             reads. Options -donotfixmates and -clear also apply.\n");
  fprintf(stdout, "  -withlibname               For -dumpfasta, -dumpnewbler and -dumpfastq, embed the libraryname in the\n");
  fprintf(stdout, "                             created files, e.g., prefix.libname.1.fastq for fastq files.\n");

  fprintf(stdout, "\n");
  if (longhelp == 0) {
    fprintf(stdout, "Use '-h' to get a discussion of what gatekeeper is.\n");
    fprintf(stdout, "\n");
  } else {
    fprintf(stdout, "The Gatekeeper ensures that data entering the assembly system meets\n");
    fprintf(stdout, "the data specification (see GateKeeper design document).  It is also\n");
    fprintf(stdout, "used for examining and partitioning the assembler data store.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Each input message is checked for semantic consistency as described in\n");
    fprintf(stdout, "the defining document for that stage.  Messages containing a UID are\n");
    fprintf(stdout, "converted to a UID,IID pair -- the assembler modules require\n");
    fprintf(stdout, "consecutive IID beginning at 1 for efficient indexing of internal and\n");
    fprintf(stdout, "disk-based data structures; gatekeeper performs this task.  Finally,\n");
    fprintf(stdout, "each message is inserted into the assembly data store.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "The GateKeeper succeeds if it consumes its entire input with less than\n");
    fprintf(stdout, "a specified number of errors (the -e option).  Upon successful exit,\n");
    fprintf(stdout, "the store reflects all of the records that were successfully read.\n");
    fprintf(stdout, "Unsuccessful records are reported to stderr, along with a brief\n");
    fprintf(stdout, "explanation of the problem.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "If unsuccessful, the store is partially updated.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Resoure Requirements\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "The key gatekeeper data structures are in-memory copies of its store.\n");
    fprintf(stdout, "This store should scale linearly with the number of fragments.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "No formal benchmarking of the gatekeeper has been performed to date.\n");
    fprintf(stdout, "However, each LKG message requires four random disk accesses -- two to\n");
    fprintf(stdout, "read the linked fragment records, and two two write the updated\n");
    fprintf(stdout, "fragment records.  This can cause problems when gatekeeper is run over\n");
    fprintf(stdout, "low-performance or heavily used NFS mount points.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "3. General Pre Conditions\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "a) Each input UID must be unique.  A new message with a duplicate UID\n");
    fprintf(stdout, "will be rejected as an error.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "b) Any object referred to in a message must be defined:\n");
    fprintf(stdout, "def-before-ref.  Likewise, before an object can be deleted, all\n");
    fprintf(stdout, "references to it must be removed: unref-before-undef.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "c) The input specification is defined elsewhere.\n");
    fprintf(stdout, "\n");
  }
}


char *
constructIIDdumpFromIDFile(char *gkpStoreName, char *iidToDump, char *uidFileName, char *iidFileName) {

  if ((uidFileName == NULL) && (iidFileName == NULL))
    return(iidToDump);

  gkStore *gkp      = new gkStore(gkpStoreName, FALSE, FALSE);

  AS_IID           lastElem = gkp->gkStore_getNumFragments() + 1;
  FILE            *F        = NULL;
  char             L[1024];

  if (iidToDump == NULL)
    iidToDump = (char *)safe_calloc(lastElem, sizeof(char));

  if (iidFileName) {
    errno = 0;
    if (strcmp(iidFileName, "-") == 0)
      F = stdin;
    else
      F = fopen(iidFileName, "r");
    if (errno) {
      fprintf(stderr, "%s: Couldn't open -iid file '%s': %s\n", progName, iidFileName, strerror(errno));
      exit(1);
    }
    fgets(L, 1024, F);
    while (!feof(F)) {
      AS_IID      iid = AS_IID_fromString(L, NULL);
      if (iid >= lastElem)
        fprintf(stderr, "%s: IID "F_IID" too big, ignored.\n", progName, iid);
      else
        iidToDump[iid]++;
      fgets(L, 1024, F);
    }
    if (F != stdin)
      fclose(F);
  }

  if (uidFileName) {
    errno = 0;
    if (strcmp(uidFileName, "-") == 0)
      F = stdin;
    else
      F = fopen(uidFileName, "r");
    if (errno) {
      fprintf(stderr, "%s: Couldn't open -uid file '%s': %s\n", progName, uidFileName, strerror(errno));
      exit(1);
    }
    fgets(L, 1024, F);
    while (!feof(F)) {
      chomp(L);
      AS_UID  uid = AS_UID_lookup(L, NULL);
      AS_IID  iid = gkp->gkStore_getUIDtoIID(uid, NULL);

      if (iid == 0)
        fprintf(stderr, "%s: UID %s doesn't exist, ignored.\n", progName, L);
      else if (iid >= lastElem)
        fprintf(stderr, "%s: UID %s is IID "F_IID", and that's too big, ignored.\n", progName, L, iid);
      else
        iidToDump[iid]++;

      fgets(L, 1024, F);
    }
    if (F != stdin)
      fclose(F);
  }

  delete gkp;

  return(iidToDump);
}



char *
constructIIDdump(char  *gkpStoreName,
                 char  *iidToDump,
                 uint32 dumpRandLib,
                 uint32 dumpRandMateNum,
                 uint32 dumpRandSingNum,
                 double dumpRandFraction,
                 uint64 dumpRandLength,
                 int	dumpInvert) {

  if ((dumpRandMateNum == 0) && (dumpRandSingNum == 0) &&
      (dumpRandFraction == 0.0) &&
      (dumpRandLength == 0))
    return(iidToDump);

  gkStore *gkp = new gkStore(gkpStoreName, FALSE, FALSE);

  uint32      numMated  = 0;
  uint32      numSingle = 0;
  uint32      numNoLib  = 0;

  uint64      lenFrag   = 0;

  uint32      numTotal  = gkp->gkStore_getNumFragments() + 1;  //  Should count, or remember this when building

  uint32     *candidatesS    = (uint32 *)safe_malloc(numTotal * sizeof(uint32));
  uint32      candidatesSLen = 0;

  uint32     *candidatesA    = (uint32 *)safe_malloc(numTotal * sizeof(uint32));
  uint32     *candidatesB    = (uint32 *)safe_malloc(numTotal * sizeof(uint32));
  uint32      candidatesMLen = 0;

  gkFragment  fr;
  gkStream   *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  if (iidToDump == NULL)
    iidToDump = (char *)safe_calloc(numTotal, sizeof(char));

  //  Scan the whole fragstore, looking for mated reads in the
  //  correct library, and save the lesser of the two reads.
  //
  while (fs->next(&fr)) {

    if (fr.gkFragment_getLibraryIID() == 0)
      numNoLib++;

    if ((dumpRandLib == 0) ||
        (dumpRandLib == fr.gkFragment_getLibraryIID())) {

		// default to selecting everything
    	if (dumpInvert) {
    		iidToDump[fr.gkFragment_getReadIID()] = 1;
    	}

      //  Build lists of singletons and mated frags in this library.
      //  Save only the smaller mate ID.

      lenFrag += fr.gkFragment_getSequenceLength();

      if (fr.gkFragment_getMateIID() == 0) {
        numSingle++;
        candidatesS[candidatesSLen] = fr.gkFragment_getReadIID();
        candidatesSLen++;
      } else if (fr.gkFragment_getReadIID() < fr.gkFragment_getMateIID()) {
        numMated += 2;
        candidatesA[candidatesMLen] = fr.gkFragment_getReadIID();
        candidatesB[candidatesMLen] = fr.gkFragment_getMateIID();
        candidatesMLen++;
      }
    }
  }

  delete fs;
  delete gkp;

  if (numNoLib)
    fprintf(stderr, "WARNING: found "F_U32" reads with no library (usually caused by using frg format 1).\n", numNoLib);

  //  Now pick N reads from our list of candidates, and let the dump
  //  routines fill in the missing mates

  srand48(srandSeed == 0 ? time(NULL) : srandSeed);

  if (dumpRandLength > 0) {
    double a = (double)lenFrag / (numSingle + numMated);
    dumpRandSingNum = (uint32)(dumpRandLength / a * numSingle / (numSingle + numMated));
    dumpRandMateNum = (uint32)(dumpRandLength / a * numMated  / (numSingle + numMated));
    //fprintf(stderr, "randLength %f %d %d\n", a, dumpRandSingNum, dumpRandMateNum);
  }

  if (dumpRandFraction > 0) {
    dumpRandSingNum = (uint32)(numSingle * dumpRandFraction);
    dumpRandMateNum = (uint32)(numMated  * dumpRandFraction);
    //fprintf(stderr, "randFraction %d %d\n", dumpRandSingNum, dumpRandMateNum);
  }

  for (uint32 i=0; (i < dumpRandSingNum) && (candidatesSLen > 0); i++) {
    int  x = lrand48() % candidatesSLen;
    iidToDump[candidatesS[x]] = !dumpInvert;
    candidatesSLen--;
    candidatesS[x] = candidatesS[candidatesSLen];
  }

  for (uint32 i=0; (i < dumpRandMateNum) && (candidatesMLen > 0); i += 2) {
    int  x = lrand48() % candidatesMLen;
    iidToDump[candidatesA[x]] = !dumpInvert;
    iidToDump[candidatesB[x]] = !dumpInvert;
    candidatesMLen--;
    candidatesA[x] = candidatesA[candidatesMLen];
    candidatesB[x] = candidatesB[candidatesMLen];
  }

  safe_free(candidatesS);
  safe_free(candidatesA);
  safe_free(candidatesB);

  return(iidToDump);
}


char *
constructIIDdumpLongest(char  *gkpStoreName,
                 char  *iidToDump,
                 uint32 dumpRandLib,
                 uint32 dumpLongestMin,
                 uint64 dumpLongestLength,
                 int dumpInvert) {

  if (dumpLongestLength == 0 && dumpLongestMin == 0) {
    return(iidToDump);
  }

  gkStore *gkp = new gkStore(gkpStoreName, FALSE, FALSE);

  uint32      numNoLib  = 0;

  uint64      lenFrag   = 0;

  uint32      numTotal  = gkp->gkStore_getNumFragments() + 1;  //  Should count, or remember this when building

  multimap<uint32, AS_IID> lenToIID;

  gkFragment  fr;
  gkStream   *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  if (iidToDump == NULL)
    iidToDump = (char *)safe_calloc(numTotal, sizeof(char));

  //  Scan the whole fragstore, looking for reads in the
  //  correct library, and save the lesser of the two reads.
  //
  while (fs->next(&fr)) {

    if (fr.gkFragment_getLibraryIID() == 0)
      numNoLib++;

    if ((dumpRandLib == 0) ||
        (dumpRandLib == fr.gkFragment_getLibraryIID())) {

	// default to selecting everything
    	if (dumpInvert) {
    		iidToDump[fr.gkFragment_getReadIID()] = 1;
    	}
      uint32 len = fr.gkFragment_getSequenceLength();
      if (len > dumpLongestMin) {
         if (dumpLongestLength == 0) {
            iidToDump[fr.gkFragment_getReadIID()] = !dumpInvert;
         } else {
            lenToIID.insert(pair<uint32, AS_IID>(len, fr.gkFragment_getReadIID()));
         }
      }
    }
  }
  delete fs;
  delete gkp;

  if (numNoLib)
    fprintf(stderr, "WARNING: found "F_U32" reads with no library (usually caused by using frg format 1).\n", numNoLib);

  if (dumpLongestLength != 0) {
    for (multimap<uint32, AS_IID>::reverse_iterator iter = lenToIID.rbegin(); iter != lenToIID.rend(); iter++) {
       iidToDump[iter->second] = !dumpInvert;
       lenFrag += iter->first;
       
       if (lenFrag >= dumpLongestLength) {
          fprintf(stderr, "Longest picked cutoff: "F_U32"\n", iter->first);
          break;
       } 
    }
  }

  return(iidToDump);
}



#define DUMP_NOTHING     0
#define DUMP_INFO        1
#define DUMP_LIBRARIES   2
#define DUMP_FRAGMENTS   3
#define DUMP_FASTA       4
#define DUMP_FRG         5
#define DUMP_NEWBLER     6
#define DUMP_LASTFRG     7
#define DUMP_FASTQ       8
#define DUMP_FEATURE     9

int
main(int argc, char **argv) {

  //  Options for everyone.  Everyone needs a GateKeeper!
  //
  char            *gkpStoreName    = NULL;

  //  Options for appending or creating:
  //
  int              append             = 0;
  int              outputExists       = 0;
  char            *vectorClearFile    = NULL;
  int              assembler          = AS_ASSEMBLER_GRANDE;
  int              firstFileArg       = 0;
  int              fixInsertSizes     = 0;
  int              packedLength       = 160;

  //  Options for partitioning
  //
  char            *partitionFile = NULL;

  //  Options for dumping:
  //
  AS_IID           begIID            = 0;
  AS_IID           endIID            = 2000000000;  //  I hope I never see an assembly with 2 billion IIDs!
  char            *uidFileName       = NULL;
  char            *iidFileName       = NULL;
  int              dump              = DUMP_NOTHING;
  int              dumpTabular       = 0;
  int              dumpWithSequence  = 0;
  int              dumpAllReads      = 0;
  int              dumpClear         = AS_READ_CLEAR_LATEST;
  int              dumpAllBases      = 0;
  int              doNotFixMates     = 0;
  int			   dumpInvert		 = 0;
  int              dumpFormat        = 2;
  char            *dumpPrefix        = NULL;
  int              dumpWithLibName   = 0;
  uint32           dumpRandLib       = 0;  //  0 means "from any library"
  uint32           dumpRandMateNum   = 0;
  uint32           dumpRandSingNum   = 0;  //  Not a command line option
  double           dumpRandFraction  = 0.0;
  uint64           dumpRandLength    = 0;
  uint32           dumpLongestMin    = 0;
  uint64           dumpLongestTotal  = 0;
  char            *iidToDump         = NULL;
  
  AS_IID           featureLibIID     = 0;
  char            *featureName       = NULL;

  int              dumpDoNotUseUIDs  = FALSE;

  progName = argv[0];
  gkpStore = NULL;
  errorFP  = stdout;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  int hlp = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      append = 1;
    } else if (strcmp(argv[arg], "-b") == 0) {
      begIID = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      endIID  = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-h") == 0) {
      hlp++;
      err++;
    } else if (strcmp(argv[arg], "-o") == 0) {
      gkpStoreName = argv[++arg];
    } else if (strcmp(argv[arg], "-v") == 0) {
      vectorClearFile = argv[++arg];
      firstFileArg    = 1;  // gets us around the input file sanity check, unused otherwise
    } else if (strcmp(argv[arg], "-T") == 0) {
      assembler = AS_ASSEMBLER_OBT;
    } else if (strcmp(argv[arg], "-F") == 0) {
      fixInsertSizes = 1;
    } else if (strcmp(argv[arg], "-P") == 0) {
      partitionFile = argv[++arg];

    } else if (strcmp(argv[arg], "-pl") == 0) {
      packedLength = atoi(argv[++arg]);

      //  Except for -b and -e, all the dump options are below.

    } else if (strcmp(argv[arg], "-tabular") == 0) {
      dumpTabular = 1;
    } else if (strcmp(argv[arg], "-uid") == 0) {
      uidFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-iid") == 0) {
      iidFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-randommated") == 0) {
      dumpRandLib      = atoi(argv[++arg]);
      dumpRandMateNum  = atoi(argv[++arg]) * 2;
      if (dumpRandMateNum == 0) {
        fprintf(stderr, "%s: -randommated told to dump 0 mates; exit.\n", argv[0]);
        exit(0);
      }
    } else if (strcmp(argv[arg], "-deterministic") == 0) {
       srandSeed = 1;

    } else if (strcmp(argv[arg], "-randomsubset") == 0) {
      dumpRandLib      = atoi(argv[++arg]);
      dumpRandFraction = atof(argv[++arg]);
      if (dumpRandFraction == 0.0) {
        fprintf(stderr, "%s: -randomsubset told to dump 0%%; exit.\n", argv[0]);
        exit(0);
      }
    } else if (strcmp(argv[arg], "-randomlength") == 0) {
      dumpRandLib      = atoi(argv[++arg]);
      dumpRandLength   = atol(argv[++arg]);
      if (dumpRandLength == 0) {
        fprintf(stderr, "%s: -randomlength told to dump 0 bases; exit.\n", argv[0]);
        exit(0);
      }
    
    } else if (strcmp(argv[arg], "-longestovermin") == 0) {
       dumpRandLib = atoi(argv[++arg]);
       dumpLongestMin = atoi(argv[++arg]);
    }
    else if (strcmp(argv[arg], "-longestlength") == 0) {
       dumpRandLib = atoi(argv[++arg]);
       dumpLongestTotal = atol(argv[++arg]);
       if (dumpLongestTotal == 0) {
        fprintf(stderr, "%s: -longesttotallength told to dump 0 bases; exit.\n", argv[0]);
        exit(0);
      }

    } else if (strcmp(argv[arg], "-dumpinfo") == 0) {
      dump = DUMP_INFO;
    } else if (strcmp(argv[arg], "-lastfragiid") == 0) {
      dump = DUMP_LASTFRG;
    } else if (strcmp(argv[arg], "-dumplibraries") == 0) {
      dump = DUMP_LIBRARIES;
    } else if (strcmp(argv[arg], "-dumpfragments") == 0) {
      dump = DUMP_FRAGMENTS;
    } else if (strcmp(argv[arg], "-withsequence") == 0) {
      dumpWithSequence = 1;
    } else if (strcmp(argv[arg], "-withlibname") == 0) {
      dumpWithLibName = 1;
    } else if (strcmp(argv[arg], "-clear") == 0) {
      dumpClear      = gkStore_decodeClearRegionLabel(argv[++arg]);
      if (dumpClear == AS_READ_CLEAR_ERROR) {
        fprintf(stderr, "%s: -clear %s is not a valid clear range.\n", argv[0], argv[arg]);
        exit(0);
      }
    } else if (strcmp(argv[arg], "-legacyformat") == 0) {
      dumpFormat = 1;
    } else if (strcmp(argv[arg], "-dumpfasta") == 0) {
      dump       = DUMP_FASTA;
      dumpPrefix = argv[++arg];
    } else if (strcmp(argv[arg], "-allreads") == 0) {
      dumpAllReads = 1;
    } else if (strcmp(argv[arg], "-allbases") == 0) {
      dumpAllBases = 1;
    } else if (strcmp(argv[arg], "-dumpfrg") == 0) {
      dump = DUMP_FRG;
    } else if (strcmp(argv[arg], "-dumpnewbler") == 0) {
      dump       = DUMP_NEWBLER;
      dumpPrefix = argv[++arg];
    } else if (strcmp(argv[arg], "-dumpfastq") == 0) {
      dump       = DUMP_FASTQ;
      dumpPrefix = argv[++arg];
    } else if (strcmp(argv[arg], "-isfeatureset") == 0 ) {
      dump = DUMP_FEATURE;
      featureLibIID = atoi(argv[++arg]);
      featureName = argv[++arg];
    } else if (strcmp(argv[arg], "-nouid") == 0 ) {
      dumpDoNotUseUIDs = 1;
    } else if (strcmp(argv[arg], "-donotfixmates") == 0) {
      doNotFixMates = 1;
    } else if (strcmp(argv[arg], "-invert") == 0) {
    	dumpInvert = 1;

        //  End of dump options, SECRET OPTIONS below

      //} else if (strcmp(argv[arg], "--rebuildmap") == 0) {
      //  gkStore_rebuildUIDtoIID(argv[arg+1]);
      //  exit(0);
    } else if (strcmp(argv[arg], "--revertclear") == 0) {
      //  Takes two args:  clear region name, gkpStore
      //
      revertClearRange(argv[arg+1], argv[arg+2]);
      exit(0);

    } else if ((strcmp(argv[arg], "--edit") == 0) ||
               (strcmp(argv[arg], "--testedit") == 0)) {
      //  Takes two args:  edit file, gkpStore
      //
      editStore(argv[arg+1], argv[arg+2], (strcmp(argv[arg], "--edit") == 0));
      exit(0);

    //  End of SECRET options

    } else if (strcmp(argv[arg], "--") == 0) {
      firstFileArg = arg++;
      arg = argc;
    } else if ((argv[arg][0] == '-') && (argv[arg][1] != 0)) {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    } else {
      firstFileArg = arg;

      if ((dump != DUMP_NOTHING) || (partitionFile))
        gkpStoreName = argv[arg];

      arg = argc;
    }

    arg++;
  }

  if ((err) || (gkpStoreName == NULL) || (firstFileArg == 0)) {
    usage(progName, hlp);
    exit(1);
  }

  if (vectorClearFile) {
    updateVectorClear(vectorClearFile, gkpStoreName);
    exit(0);
  }

  iidToDump = constructIIDdumpFromIDFile(gkpStoreName, iidToDump, uidFileName, iidFileName);
  iidToDump = constructIIDdump(gkpStoreName, iidToDump, dumpRandLib, dumpRandMateNum, dumpRandSingNum, dumpRandFraction, dumpRandLength, dumpInvert);
  iidToDump = constructIIDdumpLongest(gkpStoreName, iidToDump, dumpRandLib, dumpLongestMin, dumpLongestTotal, dumpInvert);

  if (dump != DUMP_NOTHING) {
    int exitVal = 0;
   
    switch (dump) {
      case DUMP_INFO:
        dumpGateKeeperInfo(gkpStoreName, dumpTabular);
        break;
      case DUMP_LIBRARIES:
        dumpGateKeeperLibraries(gkpStoreName, begIID, endIID, iidToDump, dumpTabular, dumpDoNotUseUIDs);
        break;
      case DUMP_FRAGMENTS:
        dumpGateKeeperFragments(gkpStoreName, begIID, endIID, iidToDump,
                                dumpWithSequence,
                                dumpClear,
                                dumpTabular,
                                dumpDoNotUseUIDs);
        break;
      case DUMP_FASTA:
        dumpGateKeeperAsFasta(gkpStoreName, dumpPrefix, dumpWithLibName, begIID, endIID, iidToDump,
                              doNotFixMates,
                              dumpAllReads, dumpAllBases,
                              dumpClear,
                              dumpDoNotUseUIDs);
        break;
      case DUMP_FRG:
        dumpGateKeeperAsFRG(gkpStoreName, dumpFormat, begIID, endIID, iidToDump,
                            doNotFixMates,
                            dumpAllReads,
                            dumpClear,
                            dumpDoNotUseUIDs);
        break;
      case DUMP_NEWBLER:
        dumpGateKeeperAsNewbler(gkpStoreName, dumpPrefix, dumpWithLibName, begIID, endIID, iidToDump,
                                doNotFixMates,
                                dumpAllReads, dumpAllBases,
                                dumpClear,
                                dumpDoNotUseUIDs);
        break;
      case DUMP_LASTFRG:
        {
          gkStore *gkp = new gkStore(gkpStoreName, FALSE, FALSE);
          fprintf(stdout, "Last frag in store is iid = %d\n", gkp->gkStore_getNumFragments());
          delete gkp;
        }
        break;
      case DUMP_FASTQ:
        dumpGateKeeperAsFastQ(gkpStoreName, dumpPrefix, dumpWithLibName, begIID, endIID, iidToDump,
                              doNotFixMates,
                              dumpAllReads, dumpAllBases,
                              dumpClear,
                              dumpDoNotUseUIDs);
        break;
      case DUMP_FEATURE:
         exitVal = (dumpGateKeeperIsFeatureSet(gkpStoreName, featureLibIID, featureName) == 0);
      default:
        break;
    }

    exit(exitVal);
  }


  if (partitionFile) {
    Build_Partition(gkpStoreName, partitionFile, GKFRAGMENT_QLT);
    exit(0);
  }


  if (append)
    //  used for updating distances after cgw
    gkpStore = new gkStore(gkpStoreName, FALSE, TRUE);
  else
    gkpStore = new gkStore(gkpStoreName, packedLength);

  //  This is a special case for gatekeeper; we never call
  //  gkStore_getFragment() and so we never set up the gkFragment.
  //  This also enables some of the set() methods.
  //
  gkFrag1 = new gkFragment();
  gkFrag2 = new gkFragment();

  gkFrag1->gkFragment_enableGatekeeperMode(gkpStore);
  gkFrag2->gkFragment_enableGatekeeperMode(gkpStore);


  sprintf(fastqUIDmapName, "%s.errorLog", gkpStoreName);

  errno = 0;
  errorFP = fopen(fastqUIDmapName, "w");
  if (errno) {
    fprintf(stderr, "%s: cannot open error file '%s': %s\n", progName, fastqUIDmapName, strerror(errno));
    exit(1);
  }

  sprintf(fastqUIDmapName, "%s.fastqUIDmap", gkpStoreName);


  for (; firstFileArg < argc; firstFileArg++) {
    GenericMesg           *pmesg             = NULL;

    fprintf(stderr, "\n");
    fprintf(stderr, "Starting file '%s'.\n", argv[firstFileArg]);

    AS_MSG_resetProtoLineNum();
    AS_MSG_setFormatVersion(1);

    compressedFileReader *inFile = new compressedFileReader(argv[firstFileArg]);

    while (EOF != ReadProtoMesg_AS(inFile->file(), &pmesg)) {
      if        (pmesg->t == MESG_DST) {
        Check_DistanceMesg((DistanceMesg *)pmesg->m, fixInsertSizes);
      } else if (pmesg->t == MESG_LIB) {
        Check_LibraryMesg((LibraryMesg *)pmesg->m, fixInsertSizes, packedLength);
      } else if (pmesg->t == MESG_FRG) {
        Check_FragMesg((FragMesg *)pmesg->m, assembler);
      } else if (pmesg->t == MESG_LKG) {
        Check_LinkMesg((LinkMesg *)pmesg->m);
      } else if (pmesg->t == MESG_PLC) {
        Check_PlacementMesg((PlacementMesg *) pmesg->m);      
      } else if (pmesg->t == MESG_VER) {
        //  Ignore
      } else {
        //  Ignore messages we don't understand
        AS_GKP_reportError(AS_GKP_UNKNOWN_MESSAGE, 0, MessageTypeName[pmesg->t]);
      }
    }

    delete inFile;
  }

  delete gkFrag1;
  delete gkFrag2;

  delete gkpStore;

  if (errorFP != stdout)
    //  Close it if it is a file
    fclose(errorFP);
  else
    //  Or flush it; it is stdout and summarizeErrors writes to stderr
    fflush(errorFP);

  if (fastqUIDmap)
    fclose(fastqUIDmap);

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");

  return(AS_GKP_summarizeErrors(gkpStoreName));
}
