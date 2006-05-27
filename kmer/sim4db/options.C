#include "sim4db.H"
#include "configuration.H"

#include "buildinfo-sim4th.h"
#include "buildinfo-libbio.h"
#include "buildinfo-libutil.h"


const char *usage = 
"usage: %s [options]\n"
"\n"
"       -v            print status to stderr while running\n"
"       -V            print script lines (stderr) as they are processed\n"
"       -YN           print script lines (to given file) as they are processed, annotated with yes/no\n"
"\n"
"       -cdna         use these cDNA sequences\n"
"       -genomic      use these genomic sequences\n"
"       -script       use this script file\n"
"       -output       write output to this file\n"
"       -stats        write execution statistics to this file\n"
"       -touch        create this file when the program finishes execution\n"
"\n"
"       -mincoverage  iteratively find all exon models with the specified\n"
"                     minimum PERCENT COVERAGE\n"
"       -minidentity  iteratively find all exon models with the specified\n"
"                     minimum PERCENT EXON IDENTITY\n"
"       -minlength    iteratively find all exon models with the specified\n"
"                     minimum ABSOLUTE COVERAGE (number of bp matched)\n"
"       -alwaysreport always report <number> exon models, even if they\n"
"                     are below the quality thresholds\n"
"\n"
"         If no mincoverage or minidentity or minlength is given, only\n"
"         the best exon model is returned.\n"
"\n"
"         You will probably want to specify ALL THREE of mincoverage,\n"
"         minidentity and minlength!  Don't assume the default values\n"
"         are what you want!\n"
"\n"
"         You will DEFINITELY want to specify at least one of mincoverage,\n"
"         minidentity and minlength with alwaysreport!  If you don't, mincoverage\n"
"         will be set to 90 and minidentity to 95 -- to reduce the number of\n"
"         spurious matches when a good match is found.\n"
"\n"
"       -nodeflines   don't include the defline in the output\n"
"       -alignments   print alignments\n"
"\n"
"       -polytails    DON'T mask poly-A and poly-T tails.\n"
"       -cut          Trim marginal exons if A/T %% > x (poly-AT tails)\n"
"\n"
"       -noncanonical Don't force canonical splice sites\n"
"\n"
"       -forcestrand  Force the strand prediction to always be\n"
"                     'forward' or 'reverse'\n"
"\n"
"       -interspecies Configure sim4 for better inter-species alignments\n"
"\n"
"  The following are for use only by immortals.\n"
"       -H            set the relink weight factor\n"
"       -K            set the first MSP threshold\n"
"       -C            set the second MSP threshold\n"
"       -Ma           set the limit of the number of MSPs allowed\n"
"       -Mp           same, as percentage of bases in cDNA\n"
"                     NOTE:  If used, both -Ma and -Mp must be specified!\n"
;




void
configuration::parseCommandLine(int argc, char **argv) {
  int arg = 1;

  while (arg < argc) {
    if        (strncmp(argv[arg], "--buildinfo", 3) == 0) {
        buildinfo_sim4th(stderr);
        buildinfo_libbio(stderr);
        buildinfo_libutil(stderr);
        exit(1);
    } else if (strncmp(argv[arg], "-alignments", 4) == 0) {
      sim4params.setPrintAlignments(true);
    } else if (strncmp(argv[arg], "-alwaysprint", 4) == 0) {
      sim4params.setFindAllExons(true);
      arg++;
      sim4params.setAlwaysReport(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-cdna", 3) == 0) {
      arg++;
      cdnaFileName = argv[arg];
    } else if (strncmp(argv[arg], "-cut", 3) == 0) {
      arg++;
      double x = atof(argv[arg]);
      if (x < 0.0) {
        fprintf(stderr, "WARNING:  -cut adjusted to 0.0 (you gave %f)!\n", x);
        x = 0.0;
      }
      if (x > 1.0) {
        fprintf(stderr, "WARNING:  -cut adjusted to 1.0 (you gave %f)!\n", x);
        x = 1.0;
      }
      sim4params.setPolyTailPercent(x);
    } else if (strncmp(argv[arg], "-genomic", 2) == 0) {
      arg++;
      databaseFileName = argv[arg];
    } else if (strncmp(argv[arg], "-minc", 5) == 0) {
      sim4params.setFindAllExons(true);
      arg++;
      sim4params.setMinCoverage(atoi(argv[arg]) / 100.0);
    } else if (strncmp(argv[arg], "-mini", 5) == 0) {
      sim4params.setFindAllExons(true);
      arg++;
      sim4params.setMinPercentExonIdentity(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-minl", 5) == 0) {
      sim4params.setFindAllExons(true);
      arg++;
      sim4params.setMinCoverageLength(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-nod", 4) == 0) {
      sim4params.setIncludeDefLine(false);
    } else if (strncmp(argv[arg], "-non", 4) == 0) {
      sim4params.setDontForceCanonicalSplicing(true);
    } else if (strncmp(argv[arg], "-f", 2) == 0) {
      sim4params.setForceStrandPrediction(true);
    } else if (strncmp(argv[arg], "-o", 2) == 0) {
      arg++;
      outputFileName = argv[arg];
    } else if (strncmp(argv[arg], "-po", 3) == 0) {
      sim4params.setIgnorePolyTails(false);
    } else if (strncmp(argv[arg], "-sc", 3) == 0) {
      arg++;
      scriptFileName = argv[arg];
    } else if (strncmp(argv[arg], "-pa", 3) == 0) {
      pairwise = true;
    } else if (strncmp(argv[arg], "-st", 3) == 0) {
      arg++;
      statsFileName = argv[arg];
    } else if (strncmp(argv[arg], "-to", 3) == 0) {
      arg++;
      touchFileName = argv[arg];
    } else if (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
    } else if (strncmp(argv[arg], "-YN", 3) == 0) {
      arg++;
      yesnoFileName = argv[arg];
    } else if (strncmp(argv[arg], "-H", 2) == 0) {
      arg++;
      sim4params.setRelinkWeight(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-K", 2) == 0) {
      arg++;
      sim4params.setMSPThreshold1(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-C", 2) == 0) {
      arg++;
      sim4params.setMSPThreshold2(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-Ma", 3) == 0) {
      arg++;
      sim4params.setMSPLimitAbsolute(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-Mp", 3) == 0) {
      arg++;
      sim4params.setMSPLimitPercent(atof(argv[arg]));
    } else if (strncmp(argv[arg], "-interspecies", 2) == 0) {
      sim4params.setInterspecies(true);
    } else {
      fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
    }

    arg++;
  }

  if (cdnaFileName == 0L) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "--No EST file?\n");
    exit(1);
  }

  if (databaseFileName == 0L) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "--No database file?\n");
    exit(1);
  }

  if (outputFileName == 0L) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "--No output file?\n");
    exit(1);
  }
}
