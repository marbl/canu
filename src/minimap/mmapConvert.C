
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Sergey Koren beginning on 2016-FEB-24
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "ovStore.H"
#include "splitToWords.H"

#include <vector>

using namespace std;


int
main(int argc, char **argv) {
  bool            asCoords = true;

  char           *outName     = NULL;

  vector<char *>  files;


  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      outName = argv[++arg];

    } else if (AS_UTL_fileExists(argv[arg])) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((err) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] file.mhap[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  Converts mhap native output to ovb\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o out.ovb     output file\n");
    fprintf(stderr, "\n");

    if (files.size() == 0)
      fprintf(stderr, "ERROR:  no overlap files supplied\n");

    exit(1);
  }

  char        *ovStr = new char [1024];

  ovOverlap   ov(NULL);
  ovFile      *of = new ovFile(outName, ovFileFullWrite);

  for (uint32 ff=0; ff<files.size(); ff++) {
    compressedFileReader  *in = new compressedFileReader(files[ff]);

    //  $1    							$2   	$3	$4 	$5  	$6 							$7   	$8   	$9  	$10 			$11  	$12	$13
    //  0     							1    	2       3   	4   	5   							6    	7    	8   	9   			10   	11	12
    //  0f1bd7b6-a7f2-4bcb-8575-d617f1394b8a_Basecall_2D_2d	8189	1310	8014	+	b74d9367-f45a-4684-8bfc-ff533629b030_Basecall_2D_2d	14205	7340	14051	277			6711	255	cm:i:32
    //  0f1bd7b6-a7f2-4bcb-8575-d617f1394b8a_Basecall_2D_2d	8189	1152	7272	-	a3026aca-57a7-4639-96bf-b76624cf2d34_Basecall_2D_2d	7731	1642	7547	157			6120	255	cm:i:24
    //  aiid  							alen    bgn	end	bori	biid 							blen	bgn	end	#match minimizers	alnlen	?	cm:i:errori
    //

    while (fgets(ovStr, 1024, in->file()) != NULL) {
      splitToWords  W(ovStr);

      ov.a_iid = W(0);
      ov.b_iid = W(5);

      if (ov.a_iid == ov.b_iid)
        continue;

      ov.dat.ovl.forUTG = true;
      ov.dat.ovl.forOBT = true;
      ov.dat.ovl.forDUP = true;

      ov.dat.ovl.ahg5 = W(2);
      ov.dat.ovl.ahg3 = W(1) - W(3);

      if (W[4][0] == '+') {
        ov.dat.ovl.bhg5 = W(7);
        ov.dat.ovl.bhg3 = W(6) - W(8);
        ov.flipped(false);
      } else {
        ov.dat.ovl.bhg3 = W(7);
        ov.dat.ovl.bhg5 = W(6) - W(8);
        ov.flipped(true);
      }

      ov.erate(1-((double)W(9)/W(10)));

      of->writeOverlap(&ov);
    }

    arg++;
  }

  delete    of;
  delete [] ovStr;

  exit(0);
}
