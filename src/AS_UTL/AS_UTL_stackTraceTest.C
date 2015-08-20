
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
 *    Brian P. Walenz on 2013-OCT-18
 *      are Copyright 2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

int crash()
{
  char * p = NULL;
  *p = 0;
  return 0;
}

int foo4()
{
  crash();
  return 0;
}

int foo3()
{
  foo4();
  return 0;
}

int foo2()
{
  foo3();
  return 0;
}

int foo1()
{
  foo2();
  return 0;
}

int main(int argc, char ** argv)
{
  struct sigaction sigact;

  sigact.sa_sigaction = crit_err_hdlr;
  sigact.sa_flags     = SA_RESTART | SA_SIGINFO;

  //  Don't especially care if this fails or not.
  sigaction(SIGSEGV, &sigact, NULL);

  foo1();

  exit(EXIT_SUCCESS);
}
