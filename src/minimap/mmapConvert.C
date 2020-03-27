
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

#include "runtime.H"
#include "ovStore.H"
#include "strings.H"

#include <vector>

using namespace std;

int
main(int argc, char **argv) {
  char           *outName  = NULL;
  char           *seqName  = NULL;
  bool		  partialOverlaps = false;
  uint32          minOverlapLength = 0;
  double          erate = 0;

  vector<char *>  files;

  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      outName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-partial") == 0) {
      partialOverlaps = true;

   } else if (strcmp(argv[arg], "-e") == 0) {
      erate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-len") == 0) {
      minOverlapLength = atoi(argv[++arg]);

    } else if (fileExists(argv[arg])) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((err) || (seqName == NULL) || (outName == NULL) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] file.mhap[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  Converts mhap native output to ovb\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o out.ovb     output file\n");
    fprintf(stderr, "\n");

    if (seqName == NULL)
      fprintf(stderr, "ERROR:  no seqStore (-S) supplied\n");
    if (files.size() == 0)
      fprintf(stderr, "ERROR:  no overlap files supplied\n");

    exit(1);
  }

  char        *ovStr = new char [1024*1024];

  sqStore    *seqStore = new sqStore(seqName);
  ovOverlap   ov;
  ovFile      *of = new ovFile(seqStore, outName, ovFileFullWrite);

  for (uint32 ff=0; ff<files.size(); ff++) {
    compressedFileReader  *in = new compressedFileReader(files[ff]);

    //  $1        $2     $3     $4     $5     $6         $7      $8    $9     $10      $11          $12        $13
    //  0         1      2      3      4      5          6       7     8      9        10           11         12
    //  aiid      alen   bgn    end    bori   biid       blen    bgn   end    #match   minimizers   alnlen     cm:i:errori
    //  read1	5064	0	5060	+	read164	7384	138	5251	4763	5144	0	tp:A:S	cm:i:1410	s1:i:4754	dv:f:0.0142
    //

    while (fgets(ovStr, 1024*1024, in->file()) != NULL) {
      splitToWords  W(ovStr);

      ov.a_iid = atoi(W[0]+4);
      ov.b_iid = atoi(W[5]+4);

      if (ov.a_iid == ov.b_iid)
        continue;

      ov.dat.ovl.ahg5 = W.toint32(2);
      ov.dat.ovl.ahg3 = W.toint32(1) - W.toint32(3);

      if (W[4][0] == '+') {
        ov.dat.ovl.bhg5 = W.toint32(7);
        ov.dat.ovl.bhg3 = W.toint32(6) - W.toint32(8);
        ov.flipped(false);
      } else {
        ov.dat.ovl.bhg3 = W.toint32(7);
        ov.dat.ovl.bhg5 = W.toint32(6) - W.toint32(8);
        ov.flipped(true);
      }

      ov.erate((double)atof(W[15]+5));

      //  Check the overlap - the hangs must be less than the read length.

      uint32  alen = seqStore->sqStore_getReadLength(ov.a_iid);
      uint32  blen = seqStore->sqStore_getReadLength(ov.b_iid);

      if ((alen < ov.dat.ovl.ahg5 + ov.dat.ovl.ahg3) ||
          (blen < ov.dat.ovl.bhg5 + ov.dat.ovl.bhg3))
        fprintf(stderr, "INVALID OVERLAP " F_U32 " (len %6d) " F_U32 " (len %6d) hangs " F_OV " " F_OV " - " F_OV " " F_OV "%s\n",
                ov.a_iid, alen,
                ov.b_iid, blen,
                ov.dat.ovl.ahg5, ov.dat.ovl.ahg3,
                ov.dat.ovl.bhg5, ov.dat.ovl.bhg3,
                (ov.dat.ovl.flipped) ? " flipped" : ""), exit(1);

      ov.dat.ovl.forUTG = (partialOverlaps == false) && (ov.overlapIsDovetail() == true);;
      ov.dat.ovl.forOBT = partialOverlaps;
      ov.dat.ovl.forDUP = partialOverlaps;

      // check the length is big enough
      if (ov.a_end() - ov.a_bgn() < minOverlapLength || ov.b_end() - ov.b_bgn() < minOverlapLength) {
         continue;
      }
      // check if the erate is OK
      if (ov.erate() > erate) {
         continue;
      }
      //  Overlap looks good, write it!

      of->writeOverlap(&ov);
    }

    arg++;
  }

  delete    of;
  delete [] ovStr;

  delete seqStore;

  exit(0);
}
