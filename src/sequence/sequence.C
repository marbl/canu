
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sequence/sequence.H"

#include "files.H"
#include "strings.H"


enum opMode {
  modeSummarize,     //  Summarize sequences in FASTA or FASTQ inputs
  modeExtract,       //  Extract sequences or subsequences from FASTA or FASTQ inputs
  modeGenerate,      //  Generate random sequences
  modeSimulate,      //  Simulate reads from FASTA or FASTQ inputs (ontigs/scaffolds/chromosomes)
  modeSample,        //  Extract random sequences from FASTA or FASTQ inputs
  modeShift,         //  Generate sequence based on a shift register
  modeMutate,        //  Randomly mutate bases
  modeUnset          //  Cause an error
};



int
main(int argc, char **argv) {
  vector<char *>              inputs;

  opMode                      mode = modeUnset;

  summarizeParameters         sumPar;
  extractParameters           extPar;
  generateParameters          genPar;
  simulateParameters          simPar;
  sampleParameters            samPar;
  shiftRegisterParameters     srPar;
  mutateParameters            mutPar;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {

    //  SUMMARIZE

    if      (strcmp(argv[arg], "summarize") == 0) {
      mode = modeSummarize;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-size") == 0)) {
      sumPar.genomeSize = strtoull(argv[++arg], NULL, 10);
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-1x") == 0)) {
      sumPar.limitTo1x = true;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-split-n") == 0)) {
      sumPar.breakAtN = true;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-simple") == 0)) {
      sumPar.asSimple = true;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-lengths") == 0)) {
      sumPar.asLength = true;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-assequences") == 0)) {
      sumPar.asSequences = true;
      sumPar.asBases     = false;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-asbases") == 0)) {
      sumPar.asSequences = false;
      sumPar.asBases     = true;
    }

    //  EXTRACT

    else if (strcmp(argv[arg], "extract") == 0) {
      mode = modeExtract;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-bases") == 0)) {
      decodeRange(argv[++arg], extPar.baseBgn, extPar.baseEnd);
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-sequences") == 0)) {
      decodeRange(argv[++arg], extPar.seqsBgn, extPar.seqsEnd);
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-reverse") == 0)) {
      extPar.asReverse = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-complement") == 0)) {
      extPar.asComplement = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-rc") == 0)) {
      extPar.asReverse = true;
      extPar.asComplement = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-upper") == 0)) {
      extPar.asUpperCase = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-lower") == 0)) {
      extPar.asLowerCase = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-length") == 0)) {
      decodeRange(argv[++arg], extPar.lensBgn, extPar.lensEnd);
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-lowermask") == 0)) {
      extPar.doMasking = true;
      extPar.maskWithN = false;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-nmask") == 0)) {
      extPar.doMasking = true;
      extPar.maskWithN = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "") == 0)) {
    }


    //  GENERATE

    else if (strcmp(argv[arg], "generate") == 0) {
      mode = modeGenerate;
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-min") == 0)) {
      genPar.minLength = strtouint64(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-max") == 0)) {
      genPar.maxLength = strtouint64(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-sequences") == 0)) {
      genPar.nSeqs = strtouint64(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-bases") == 0)) {
      genPar.nBases = strtouint64(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-guassian") == 0)) {
      genPar.useGaussian = true;
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-mirror") == 0)) {
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-gc") == 0)) {
      double  gc = strtodouble(argv[++arg]);
      double  at = 1.0 - gc;

      genPar.gFreq = genPar.cFreq = gc / 2.0;
      genPar.aFreq = genPar.tFreq = at / 2.0;
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-at") == 0)) {
      double  at = strtodouble(argv[++arg]);
      double  gc = 1.0 - at;

      genPar.gFreq = genPar.cFreq = gc / 2.0;
      genPar.aFreq = genPar.tFreq = at / 2.0;
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-a") == 0) ){
      genPar.aFreq = strtodouble(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-c") == 0)) {
      genPar.cFreq = strtodouble(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-g") == 0)) {
      genPar.gFreq = strtodouble(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-t") == 0)) {
      genPar.tFreq = strtodouble(argv[++arg]);
    }

    //  SIMULATE

    else if (strcmp(argv[arg], "simulate") == 0) {
      mode = modeSimulate;
    }

    else if ((mode == modeSimulate) && (strcmp(argv[arg], "-genomesize") == 0)) {
      simPar.genomeSize = strtouint64(argv[++arg]);
    }

    else if ((mode == modeSimulate) && (strcmp(argv[arg], "-coverage") == 0)) {
      simPar.desiredCoverage = strtodouble(argv[++arg]);
    }

    else if ((mode == modeSimulate) && (strcmp(argv[arg], "-nreads") == 0)) {
      simPar.desiredNumReads = strtouint64(argv[++arg]);
    }

    else if ((mode == modeSimulate) && (strcmp(argv[arg], "-nbases") == 0)) {
      simPar.desiredNumBases = strtouint64(argv[++arg]);
    }

    else if ((mode == modeSimulate) && (strcmp(argv[arg], "-circular") == 0)) {
      simPar.circular = true;
    }

    else if ((mode == modeSimulate) && (strcmp(argv[arg], "-genome") == 0)) {
      strncpy(simPar.genomeName, argv[++arg], FILENAME_MAX);
    }

    else if ((mode == modeSimulate) && (strcmp(argv[arg], "-distribution") == 0)) {
      strncpy(simPar.distribName, argv[++arg], FILENAME_MAX);
    }

    else if ((mode == modeSimulate) && (strcmp(argv[arg], "-length") == 0)) {
      decodeRange(argv[++arg], simPar.desiredMinLength, simPar.desiredMaxLength);
    }

    //else if ((mode == modeSimulate) && (strcmp(argv[arg], "-name") == 0)) {
    //  strncpy(simPar.sequenceName, argv[++arg], FILENAME_MAX);
    //}

    else if ((mode == modeSimulate) && (strcmp(argv[arg], "-output") == 0)) {
      strncpy(simPar.outputName, argv[++arg], FILENAME_MAX);
    }


    //  SAMPLE

    else if (strcmp(argv[arg], "sample") == 0) {
      mode = modeSample;
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-paired") == 0)) {
      samPar.isPaired = true;
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-output") == 0)) {
      strncpy(samPar.output1, argv[++arg], FILENAME_MAX);  //  #'s in the name will be replaced
      strncpy(samPar.output2, argv[  arg], FILENAME_MAX);  //  by '1' or '2' later.
    }


    else if ((mode == modeSample) && (strcmp(argv[arg], "-coverage") == 0)) {      //  Sample reads up to some coverage C
      samPar.desiredCoverage = strtodouble(argv[++arg]);
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-genomesize") == 0)) {
      samPar.genomeSize = strtouint64(argv[++arg]);
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-bases") == 0)) {         //  Sample B bases
      samPar.desiredNumBases = strtouint64(argv[++arg]);
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-reads") == 0)) {         //  Sample N reads
      samPar.desiredNumReads = strtouint64(argv[++arg]);
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-pairs") == 0)) {         //  Sample N pairs of reads
      samPar.desiredNumReads = strtouint64(argv[++arg]) * 2;
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-fraction") == 0)) {      //  Sample F fraction
      samPar.desiredFraction = strtodouble(argv[++arg]);
    }

    //  SHIFT

    else if (strcmp(argv[arg], "shift") == 0) {
      mode = modeShift;
    }

    else if ((mode == modeShift) && (strcmp(argv[arg], "-len") == 0)) {
    }

    else if ((mode == modeShift) && (strcmp(argv[arg], "-init") == 0)) {
      strcpy(srPar.sr, argv[++arg]);
    }

    else if ((mode == modeShift) && (strcmp(argv[arg], "-map") == 0)) {
      strcpy(srPar.sv, argv[++arg]);
    }

    else if ((mode == modeShift) && (strcmp(argv[arg], "") == 0)) {
    }

    //  MUTATE

    else if (strcmp(argv[arg], "mutate") == 0) {
      mode = modeMutate;
    }

    else if ((mode == modeMutate) && (strcmp(argv[arg], "-s") == 0)) {
      double  p = strtodouble(argv[arg+1]);
      char    a = argv[arg+2][0];
      char    b = argv[arg+3][0];

      arg += 3;

      mutPar.setProbabilitySubstititue(p, a, b);
    }

    else if ((mode == modeMutate) && (strcmp(argv[arg], "-i") == 0)) {
      double  p = strtodouble(argv[arg+1]);
      char    a = argv[arg+2][0];

      arg += 2;

      mutPar.setProbabilityInsert(p, a);
    }

    else if ((mode == modeMutate) && (strcmp(argv[arg], "-d") == 0)) {
      double  p = strtodouble(argv[arg+1]);
      char    a = argv[arg+2][0];

      arg += 2;

      mutPar.setProbabilityDelete(p, a);
    }

    //  INPUTS

    else if (fileExists(argv[arg]) == true) {
      inputs.push_back(argv[arg]);
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown parameter '%s'\n", argv[arg]);
      err.push_back(s);
    }
    arg++;
  }

  //  Check for required options.

  if (mode == modeUnset) {
    err.push_back("ERROR:  No mode (summarize, extract, generate or simulate) specified.\n");
  }
  if  (mode == modeSummarize) {
    if (inputs.size() == 0)
      err.push_back("ERROR:  No input sequence files supplied.\n");
  }
  if  (mode == modeExtract) {
    if (inputs.size() == 0)
      err.push_back("ERROR:  No input sequence files supplied.\n");
  }
  if  (mode == modeGenerate) {
  }
  if  (mode == modeSimulate) {
  }
  if  (mode == modeSample) {
  }
  if  (mode == modeShift) {
  }
  if  (mode == modeMutate) {
  }

  //  If errors, report usage.

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [mode] [options] [sequence_file ...]\n", argv[0]);
    fprintf(stderr, "\n");

    if (mode == modeUnset) {
      fprintf(stderr, "MODES:\n");
      fprintf(stderr, "  summarize      report N50, length histogram, mono-, di- and tri-nucleotide frequencies\n");
      fprintf(stderr, "  extract        extract the specified sequences\n");
      fprintf(stderr, "  sample         emit existing sequences randomly\n");
      fprintf(stderr, "  generate       generate random sequences\n");
      fprintf(stderr, "  simulate       errors in existing sequences\n");
      fprintf(stderr, "\n");
    }

    if ((mode == modeUnset) || (mode == modeSummarize)) {
      fprintf(stderr, "OPTIONS for summarize mode:\n");
      fprintf(stderr, "  -size          base size to use for N50 statistics\n");
      fprintf(stderr, "  -1x            limit NG table to 1x coverage\n");
      fprintf(stderr, "  -assequences   load data as complete sequences (for testing)\n");
      fprintf(stderr, "  -asbases       load data as blocks of bases    (for testing)\n");
      fprintf(stderr, "\n");
    }

    if ((mode == modeUnset) || (mode == modeExtract)) {
      fprintf(stderr, "OPTIONS for extract mode:\n");
      fprintf(stderr, "  -bases     baselist extract bases as specified in the 'list' from each sequence\n");
      fprintf(stderr, "  -sequences seqlist  extract ordinal sequences as specified in the 'list'\n");
      fprintf(stderr, "  -reverse            reverse the bases in the sequence\n");
      fprintf(stderr, "  -complement         complement the bases in the sequence\n");
      fprintf(stderr, "  -rc                 alias for -reverse -complement\n");
      fprintf(stderr, "  -upcase\n");
      fprintf(stderr, "  -downcase\n");
      fprintf(stderr, "  -length min-max     print sequence if it is at least 'min' bases and at most 'max' bases long\n");
      fprintf(stderr, "  \n");
      fprintf(stderr, "                      a 'baselist' is a set of integers formed from any combination\n");
      fprintf(stderr, "                      of the following, seperated by a comma:\n");
      fprintf(stderr, "                           num       a single number\n");
      fprintf(stderr, "                           bgn-end   a range of numbers:  bgn <= end\n");
      fprintf(stderr, "                      bases are spaced-based; -bases 0-2,4 will print the bases between\n");
      fprintf(stderr, "                      the first two spaces (the first two bases) and the base after the\n");
      fprintf(stderr, "                      fourth space (the fifth base).\n");
      fprintf(stderr, "  \n");
      fprintf(stderr, "                      a 'seqlist' is a set of integers formed from any combination\n");
      fprintf(stderr, "                      of the following, seperated by a comma:\n");
      fprintf(stderr, "                           num       a single number\n");
      fprintf(stderr, "                           bgn-end   a range of numbers:  bgn <= end\n");
      fprintf(stderr, "                      sequences are 1-based; -sequences 1,3-5 will print the first, third,\n");
      fprintf(stderr, "                      fourth and fifth sequences.\n");
      fprintf(stderr, "  \n");
    }

    if ((mode == modeUnset) || (mode == modeSimulate)) {
      fprintf(stderr, "OPTIONS for simulate mode:\n");
      fprintf(stderr, "  -genome G           sample reads from these sequences\n");
      fprintf(stderr, "  -circular           threat the sequences in G as circular\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -genomesize g       genome size to use for deciding coverage below\n");
      fprintf(stderr, "  -coverage c         generate approximately c coverage of output\n");
      fprintf(stderr, "  -nreads n           generate exactly n reads of output\n");
      fprintf(stderr, "  -nbases n           generate approximately n bases of output\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -distribution F     generate read length by sampling the distribution in file F\n");
      fprintf(stderr, "                        one column  - each line is the length of a sequence\n");
      fprintf(stderr, "                        two columns - each line has the 'length' and 'number of sequences'\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -length min[-max]   (not implemented)\n");
      fprintf(stderr, "  -output x.fasta     (not implemented)\n");
      fprintf(stderr, "\n");
    }

    if ((mode == modeUnset) || (mode == modeSample)) {
      fprintf(stderr, "OPTIONS for sample mode:\n");
      fprintf(stderr, "  -paired             treat inputs as paired sequences; the first two files form the\n");
      fprintf(stderr, "                      first pair, and so on.\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -output O           write output sequences to file O.  If paired, two files must be supplied.\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -coverage C         output C coverage of sequences, based on genome size G.\n");
      fprintf(stderr, "  -genomesize G       \n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -bases B            output B bases.\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -reads R            output R reads.\n");
      fprintf(stderr, "  -pairs P            output P pairs (only if -paired).\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -fraction F         output fraction F of the input bases.\n");
      fprintf(stderr, "\n");
    }

    if ((mode == modeUnset) || (mode == modeGenerate)) {
      fprintf(stderr, "OPTIONS for generate mode:\n");
      fprintf(stderr, "  -min M         minimum sequence length\n");
      fprintf(stderr, "  -max M         maximum sequence length\n");
      fprintf(stderr, "  -sequences N   generate N sequences\n");
      fprintf(stderr, "  -bases B       generate at least B bases, no more than B+maxLength-1 bases.\n");
      fprintf(stderr, "  -gaussian      99.73%% of the reads (3 standard deviations) will be between min and max\n");
      fprintf(stderr, "  -mirror F      \n");
      fprintf(stderr, "  -gc bias       sets GC/AT composition (default 0.50)\n");
      fprintf(stderr, "  -at bias       sets GC/AT composition (default 0.50)\n");
      fprintf(stderr, "  -a freq        sets frequency of A bases (default 0.25)\n");
      fprintf(stderr, "  -c freq        sets frequency of C bases (default 0.25)\n");
      fprintf(stderr, "  -g freq        sets frequency of G bases (default 0.25)\n");
      fprintf(stderr, "  -t freq        sets frequency of T bases (default 0.25)\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "The -gc option is a shortcut for setting all four base frequencies at once.  Order matters!\n");
      fprintf(stderr, "  -gc 0.6 -a 0.1 -t 0.3 -- sets G = C = 0.3, A = 0.1, T = 0.3\n");
      fprintf(stderr, "  -a 0.1 -t 0.3 -gc 0.6 -- sets G = C = 0.3, A = T = 0.15\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "Base frequencies are scaled to sum to 1.0.\n");
      fprintf(stderr, "  -a 1.25 -- results in a sum of 2.0 (1.25 + 0.25 + 0.25 + 0.25) so final frequencies will be:\n");
      fprintf(stderr, "             A =         1.25/2 = 0.625\n");
      fprintf(stderr, "             C = G = T = 0.25/2 = 0.125.\n");
      fprintf(stderr, "  -gc 0.8 -a 1.0 -t 0.2 -- sum is also 2.0, final frequencies will be:\n");
      fprintf(stderr, "             A =         1.00/2 = 0.5\n");
      fprintf(stderr, "             C = G =     0.40/2 = 0.2\n");
      fprintf(stderr, "             T =         0.20/2 = 0.1\n");
      fprintf(stderr, "\n");
    }

    if ((mode == modeUnset) || (mode == modeSimulate)) {
      fprintf(stderr, "OPTIONS for simulate mode:\n");
      fprintf(stderr, "\n");
    }

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  sumPar.finalize();
  genPar.finalize();
  simPar.finalize();
  extPar.finalize();

  switch (mode) {
    case modeSummarize:
      doSummarize(inputs, sumPar);
      break;
    case modeExtract:
      doExtract(inputs, extPar);
      break;
    case modeGenerate:
      doGenerate(genPar);
      break;
    case modeSimulate:
      doSimulate(inputs, simPar);
      break;
    case modeSample:
      doSample(inputs, samPar);
      break;
    case modeShift:
      doShiftRegister(srPar);
      break;
    case modeMutate:
      doMutate(inputs, mutPar);
      break;
    default:
      break;
  }

  return(0);
}
