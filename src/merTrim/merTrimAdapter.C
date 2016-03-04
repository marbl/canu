
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
 *  This file is derived from:
 *
 *    src/AS_MER/merTrimAdapter.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-MAY-10 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "bio++.H"
#include "sweatShop.H"
#include "existDB.H"
#include "positionDB.H"
#include "libmeryl.H"
#include "logMsg.H"

//  Size of the begin/end substring to use when creating end-to-end sequences.  This must be larger
//  than 2*merSize-1, so there is a kmer spanning the junction on both sides.  Larger is not a
//  problem; it slightly slows down the table build.
//
#define KSIZE 32


char *
createAdapterString(bool adapIllumina, bool adap454) {
  uint32  na = 0;

  uint32         adapterMax = 256;
  const char    *adapterNam[256];
  const char    *adapterSeq[256];

  if (adapIllumina) {
    adapterNam[na] = "Illumina Single End Adapter 1";                   adapterSeq[na++] = "ACACTCTTTCCCTACACGACGCTGTTCCATCT";
    adapterNam[na] = "Illumina Single End Adapter 2";                   adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Single End PCR Primer 1";                adapterSeq[na++] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Single End PCR Primer 2";                adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Single End Sequencing Primer";           adapterSeq[na++] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";

    adapterNam[na] = "Illumina Paired End Adapter 1";                   adapterSeq[na++] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Paired End Adapter 2";                   adapterSeq[na++] = "CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Paried End PCR Primer 1";                adapterSeq[na++] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Paired End PCR Primer 2";                adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Paried End Sequencing Primer 1";         adapterSeq[na++] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Paired End Sequencing Primer 2";         adapterSeq[na++] = "CGGTCTCGGCATTCCTACTGAACCGCTCTTCCGATCT";

    adapterNam[na] = "Illumina DpnII expression Adapter 1";             adapterSeq[na++] = "ACAGGTTCAGAGTTCTACAGTCCGAC";
    adapterNam[na] = "Illumina DpnII expression Adapter 2";             adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina DpnII expression PCR Primer 1";          adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina DpnII expression PCR Primer 2";          adapterSeq[na++] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA";
    adapterNam[na] = "Illumina DpnII expression Sequencing Primer";     adapterSeq[na++] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC";

    adapterNam[na] = "Illumina NlaIII expression Adapter 1";            adapterSeq[na++] = "ACAGGTTCAGAGTTCTACAGTCCGACATG";
    adapterNam[na] = "Illumina NlaIII expression Adapter 2";            adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina NlaIII expression PCR Primer 1";         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina NlaIII expression PCR Primer 2";         adapterSeq[na++] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA";
    adapterNam[na] = "Illumina NlaIII expression Sequencing Primer";    adapterSeq[na++] = "CCGACAGGTTCAGAGTTCTACAGTCCGACATG";

    adapterNam[na] = "Illumina Small RNA Adapter 1";                    adapterSeq[na++] = "GTTCAGAGTTCTACAGTCCGACGATC";
    adapterNam[na] = "Illumina Small RNA Adapter 2";                    adapterSeq[na++] = "TCGTATGCCGTCTTCTGCTTGT";
    adapterNam[na] = "Illumina Small RNA RT Primer";                    adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina Small RNA PCR Primer 1";                 adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina Small RNA PCR Primer 2";                 adapterSeq[na++] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA";
    adapterNam[na] = "Illumina Small RNA Sequencing Primer";            adapterSeq[na++] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC";

    adapterNam[na] = "Illumina Multiplexing Adapter 1";                 adapterSeq[na++] = "GATCGGAAGAGCACACGTCT";
    adapterNam[na] = "Illumina Multiplexing Adapter 2";                 adapterSeq[na++] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Multiplexing PCR Primer 1.01";           adapterSeq[na++] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Multiplexing PCR Primer 2.01";           adapterSeq[na++] = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Multiplexing Read1 Sequencing Primer";   adapterSeq[na++] = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    adapterNam[na] = "Illumina Multiplexing Index Sequencing Primer";   adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
    adapterNam[na] = "Illumina Multiplexing Read2 Sequencing Primer";   adapterSeq[na++] = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";

    adapterNam[na] = "Illumina PCR Primer Index 1";                     adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 2";                     adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 3";                     adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 4";                     adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 5";                     adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 6";                     adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 7";                     adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 8";                     adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 9";                     adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 10";                    adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 11";                    adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC";
    adapterNam[na] = "Illumina PCR Primer Index 12";                    adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC";

    adapterNam[na] = "Illumina DpnII Gex Adapter 1";                    adapterSeq[na++] = "GATCGTCGGACTGTAGAACTCTGAAC";
    adapterNam[na] = "Illumina DpnII Gex Adapter 1.01";                 adapterSeq[na++] = "ACAGGTTCAGAGTTCTACAGTCCGAC";
    adapterNam[na] = "Illumina DpnII Gex Adapter 2";                    adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina DpnII Gex Adapter 2.01";                 adapterSeq[na++] = "TCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "Illumina DpnII Gex PCR Primer 1";                 adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina DpnII Gex PCR Primer 2";                 adapterSeq[na++] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA";
    adapterNam[na] = "Illumina DpnII Gex Sequencing Primer";            adapterSeq[na++] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC";

    adapterNam[na] = "Illumina NlaIII Gex Adapter 1.01";                adapterSeq[na++] = "TCGGACTGTAGAACTCTGAAC";
    adapterNam[na] = "Illumina NlaIII Gex Adapter 1.02";                adapterSeq[na++] = "ACAGGTTCAGAGTTCTACAGTCCGACATG";
    adapterNam[na] = "Illumina NlaIII Gex Adapter 2.01";                adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina NlaIII Gex Adapter 2.02";                adapterSeq[na++] = "TCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "Illumina NlaIII Gex PCR Primer 1";                adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina NlaIII Gex PCR Primer 2";                adapterSeq[na++] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA";
    adapterNam[na] = "Illumina NlaIII Gex Sequencing Primer";           adapterSeq[na++] = "CCGACAGGTTCAGAGTTCTACAGTCCGACATG";

    adapterNam[na] = "Illumina Small RNA RT Primer";                    adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina 5p RNA Adapter";                         adapterSeq[na++] = "GTTCAGAGTTCTACAGTCCGACGATC";
    adapterNam[na] = "Illumina RNA Adapter1";                           adapterSeq[na++] = "TCGTATGCCGTCTTCTGCTTGT";

    adapterNam[na] = "Illumina Small RNA 3p Adapter 1";                 adapterSeq[na++] = "ATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "Illumina Small RNA PCR Primer 1";                 adapterSeq[na++] = "CAAGCAGAAGACGGCATACGA";
    adapterNam[na] = "Illumina Small RNA PCR Primer 2";                 adapterSeq[na++] = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA";
    adapterNam[na] = "Illumina Small RNA Sequencing Primer";            adapterSeq[na++] = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC";

    adapterNam[na] = "TruSeq Universal Adapter";                        adapterSeq[na++] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    adapterNam[na] = "TruSeq Adapter, Index 1";                         adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 2";                         adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 3";                         adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 4";                         adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 5";                         adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 6";                         adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 7";                         adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 8";                         adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 9";                         adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 10";                        adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 11";                        adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG";
    adapterNam[na] = "TruSeq Adapter, Index 12";                        adapterSeq[na++] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG";

    adapterNam[na] = "Illumina RNA RT Primer";                          adapterSeq[na++] = "GCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "Illumina RNA PCR Primer";                         adapterSeq[na++] = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA";

    adapterNam[na] = "RNA PCR Primer, Index 1";                         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 2";                         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 3";                         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 4";                         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 5";                         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 6";                         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 7";                         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 8";                         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 9";                         adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 10";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 11";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 12";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 13";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 14";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 15";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 16";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 17";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 18";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 19";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 20";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 21";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 22";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 23";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 24";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 25";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 26";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 27";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 28";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 29";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 30";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 31";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 32";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 33";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 34";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 35";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 36";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 37";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 38";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 39";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 40";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 41";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 42";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 43";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 44";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 45";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 46";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 47";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
    adapterNam[na] = "RNA PCR Primer, Index 48";                        adapterSeq[na++] = "CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";

    //  Cre-Lox from Van Nieuwerburgh, et al. "Illumina mate-paired DNA sequencing-library preparation
    //  using Cre-Lox recombination", Nucleic Acids Research, 2012, Vol 40, No. 3

    adapterNam[na] = "CreLox 5'3'";                                     adapterSeq[na++] = "CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA";
    adapterNam[na] = "CreLox 3'5'";                                     adapterSeq[na++] = "GCATTATTGAAGCATATCGTATGTAATATGCTTCAATATGCT";
  }

  if (0) {
    adapterNam[na] = "ABI Dynabead EcoP Oligo";                         adapterSeq[na++] = "CTGATCTAGAGGTACCGGATCCCAGCAGT";
    adapterNam[na] = "ABI Solid3 Adapter A";                            adapterSeq[na++] = "CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG";
    adapterNam[na] = "ABI Solid3 Adapter B";                            adapterSeq[na++] = "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT";
    adapterNam[na] = "ABI Solid3 5' AMP Primer";                        adapterSeq[na++] = "CCACTACGCCTCCGCTTTCCTCTCTATG";
    adapterNam[na] = "ABI Solid3 3' AMP Primer";                        adapterSeq[na++] = "CTGCCCCGGGTTCCTCATTCT";
    adapterNam[na] = "ABI Solid3 EF1 alpha Sense Primer";               adapterSeq[na++] = "CATGTGTGTTGAGAGCTTC";

    adapterNam[na] = "ABI Solid3 EF1 alpha Antisense Primer";           adapterSeq[na++] = "GAAAACCAAAGTGGTCCAC";
    adapterNam[na] = "ABI Solid3 GAPDH Forward Primer";                 adapterSeq[na++] = "TTAGCACCCCTGGCCAAGG";
    adapterNam[na] = "ABI Solid3 GAPDH Reverse Primer";                 adapterSeq[na++] = "CTTACTCCTTGGAGGCCATG";
  }

  if (adap454) {
    adapterNam[na] = "454 FLX Linker";                                  adapterSeq[na++] = "GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC";
    adapterNam[na] = "454 Titanium Linker";                             adapterSeq[na++] = "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG";

    adapterNam[na] = "454 AdaptorA";                                    adapterSeq[na++] = "CTGAGACAGGGAGGGAACAGATGGGACACGCAGGGATGAGATGG";
    adapterNam[na] = "454 AdaptorB";                                    adapterSeq[na++] = "CTGAGACACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGG";
  }

  assert(na < 255);

  //  Build a new sequence, from every adapter above, then every pair of adapter, F-F, F-R, R-F, R-R.
  //
  //  For the pairs, we only care about the junction, so can limit to the two kmers on either side.

  uint32   adapterListMax = 16 * 1024 * 1024;
  char    *adapterList = new char [adapterListMax];
  char    *adapterEnd  = adapterList;

  memset(adapterList, 0, sizeof(char) * 16 * 1024 * 1024);

  int32 ladapter[256];

  memset(ladapter, 0, sizeof(int32) * 256);

  char  fadapterA[256];
  char  radapterA[256];

  char  fadapterB[256];
  char  radapterB[256];


  for (uint32 a=0; a<na; a++) {
    ladapter[a] = strlen(adapterSeq[a]);

    *adapterEnd++ = '>';
    *adapterEnd++ = '\n';
    strcat(adapterEnd, adapterSeq[a]);

    while (*++adapterEnd)
      ;

    //*adapterEnd++ = 'N';
    *adapterEnd++ = '\n';
    *adapterEnd   = 0;
  }

  //  Iterate over all pairs of adapter.  For each 'a', save the last KSIZE bases.  The first/last KSIZE
  //  bases of 'b' is then appended, and this string is added to our list.

  for (uint32 a=0; a<na; a++) {
    strcpy(fadapterA, adapterSeq[a] + ((ladapter[a] < KSIZE) ? (0) : (ladapter[a] - KSIZE)));

    for (uint32 b=0; b<na; b++) {
      strcpy(fadapterB, adapterSeq[b]);

      fadapterB[KSIZE] = 0;  //  If shorter than KSIZE, already 0 terminated

      *adapterEnd++ = '>';
      *adapterEnd++ = '\n';
      *adapterEnd   = 0;
      strcat(adapterEnd, fadapterA);
      strcat(adapterEnd, fadapterB);

      while (*++adapterEnd)
        ;

      //*adapterEnd++ = 'N';
      *adapterEnd++ = '\n';
      *adapterEnd   = 0;
    }
  }

  assert(adapterEnd < adapterList + adapterListMax);

  return(adapterList);
}
