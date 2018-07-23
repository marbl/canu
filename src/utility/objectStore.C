
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "objectStore.H"
#include "strings.H"

#include <libgen.h>
#include <sys/wait.h>



extern char **environ;  //  Where, or where, is this really defined?!



static
char *
findSeqStorePath(char *requested) {
  splitToWords  F(requested, splitPaths);

  if (F.numWords() < 2)
    return(NULL);

  char  *filename  = F.last(0);
  char  *storename = F.last(1);

  //  If not a blobs file name, return no file.

  if (strlen(filename) != 10)
    return(NULL);

  if ((filename[0] != 'b') ||
      (filename[1] != 'l') ||
      (filename[2] != 'o') ||
      (filename[3] != 'b') ||
      (filename[4] != 's') ||
      (filename[5] != '.') ||
      (isdigit(filename[6]) == 0) ||
      (isdigit(filename[7]) == 0) ||
      (isdigit(filename[8]) == 0) ||
      (isdigit(filename[9]) == 0))
    return(NULL);

  //  Now just paste the two components together in the proper
  //  way and return it.

  char  *filepath = new char [FILENAME_MAX + 1];

  snprintf(filepath, FILENAME_MAX, "%s/%s", storename, filename);

  return(filepath);
}



static
char *
findOvlStorePath(char *requested) {
  splitToWords  F(requested, splitPaths);

  if (F.numWords() < 2)
    return(NULL);

  char  *basename  = NULL;
  char  *storename = F.last(1);
  char  *filename  = F.last(0);

  if (strlen(filename) != 9)
    return(NULL);

  //  If not an overlap store data file name, return no file.

  if ((isdigit(filename[0]) == 0) ||
      (isdigit(filename[1]) == 0) ||
      (isdigit(filename[2]) == 0) ||
      (isdigit(filename[3]) == 0) ||
      (filename[4]          != '<') ||
      (isdigit(filename[5]) == 0) ||
      (isdigit(filename[6]) == 0) ||
      (isdigit(filename[7]) == 0) ||
      (filename[8]          != '>'))
    return(NULL);

  //  Get ready for some ugly string parsing.  We expect strings similar to:
  //
  //    requested file F -- '../asm.ovlStore/0001<000>'
  //    current path   P -- '/path/to/assembly/correction/2-correction'
  //
  //  If the first component of F is '..', we drop it and the last component of P.
  //  When there are no more '..'s at the start, we should be left with the
  //  store name in F and the assembly stage in P.

  char  *cwd = getcwd(new char [FILENAME_MAX+1], FILENAME_MAX);

  splitToWords  P(cwd, splitPaths);

  delete [] cwd;

  uint32  nStrip = 0;

  //fprintf(stderr, "FROM cwd       '%s'\n", cwd);
  //fprintf(stderr, "     requested '%s'\n", requested);

  //  Remove identity components.

  while ((F.numWords() > 0) &&
         (strcmp(F.first(), ".") == 0))
    F.shift();

  //  Remove up components.

  while ((P.numWords() > 0) &&
         (F.numWords() > 0) &&
         (strcmp(F.first(), "..") == 0)) {
    //fprintf(stderr, "STRIP '%s' from requested and '%s' from cwd\n", F.first(), P.last());

    F.shift();
    P.pop();

    nStrip++;
  }

  //fprintf(stderr, "P.last  '%s'\n", P.last());
  //fprintf(stderr, "F.first '%s'\n", F.first());

  //  We can run in one of three different places:
  //    1) assembly_root/correction/1-stuff  - ../asm.ovlStore/0001<001>
  //    2) assembly_root/correction          -  ./asm.ovlStore/0001<001>
  //    3) assembly_root                     -  ./correction/asm.ovlStore/0001<001>
  //
  //  In the first case, we strip off the '..' and '1-stuff', set basename
  //  to the last component in P, the storename to the first component
  //  in F and the file to the last component in F (which is always true).
  //
  //  In the second case, nothing was stripped, and the result is the same.
  //
  //  In the third case, again, nothing was stripped, but the basename is
  //  now in F, not P.
  //
  //  All that boils down to

  if      (nStrip > 0) {                 //  First case.
    basename  = P.last();
    storename = F.first();
    assert(F.numWords() == 2);
  }

  else if (F.numWords() == 2) {          //  Second case.
    basename  = P.last();                //  (same result as third case)
    storename = F.first();
    assert(F.numWords() == 2);
  }

  else {                                 //  Third case.
    basename  = F.first(0);
    storename = F.first(1);
    assert(F.numWords() == 3);
  }

  //  We could check that the namespace -- the name of this assembly -- is before
  //  the basename (lots of work) and that the basename is one of 'correction',
  // 'trimming', etc.  But why?

  char  *filepath = new char [FILENAME_MAX + 1];

  //fprintf(stderr, "MAKE PATH STAGE     '%s'\n", basename);
  //fprintf(stderr, "          STORENAME '%s'\n", storename);
  //fprintf(stderr, "          FILENAME  '%s'\n", filename);

  snprintf(filepath, FILENAME_MAX, "%s/%s/%s", basename, storename, filename);

  return(filepath);
}



void
fetchFromObjectStore(char *requested) {

  //  Decide if we even need to bother.  If the file exists locally, or if
  //  one of the environment variables is missing, no, we don't need to bother.

  if (fileExists(requested))
    return;

  char  *dx = getenv("CANU_OBJECT_STORE_CLIENT");
  char  *ns = getenv("CANU_OBJECT_STORE_NAMESPACE");
  char  *pr = getenv("CANU_OBJECT_STORE_PROJECT");

  if ((dx == NULL) ||
      (ns == NULL) ||
      (pr == NULL))
    return;

  //  But now we've got a bit of a problem.  Paths to stores are relative, but we need
  //  them rooted in the assembly root directory:
  //      ../../asm.seqStore -> ./asm.seqStore
  //      ../asm.ovlStore    -> ./correction/asm.ovlStore
  //
  //  For the seqStore, we can just grab the last two components.
  //
  //  For the ovlStore, we need to parse out the subdirectory the store is in.

  fprintf(stderr, "fetchFromObjectStore()-- requesting file '%s'\n", requested);

  char *path = NULL;

  if (path == NULL)
    path = findSeqStorePath(requested);

  if (path == NULL)
    path = findOvlStorePath(requested);

  if (path == NULL)
    fprintf(stderr, "fetchFromObjectStore()-- requested file '%s', but don't know where that is.\n", requested), exit(1);

  fprintf(stderr, "fetchFromObjectStore()-- found path '%s'\n", path);

  char *cmd = new char [FILENAME_MAX+1];
  snprintf(cmd, FILENAME_MAX, "%s:%s/%s", pr, ns, path);
  char *args[8] = {"dx", "download", "--overwrite", "--no-progress", "--output", "", "", (char*)0};
  args[5] = requested;
  args[6] = cmd;

  fprintf(stderr, "fetchFromObjectStore()-- executing '%s'\n", cmd);

  int32 err;
  int32 pid = vfork();
  if ( pid == -1)
    fprintf(stderr, "vfork failed with error '%s'.\n", strerror(errno));

  if ( pid == 0 ) {
    execve(dx, args, environ);
    fprintf(stderr, "execve failed with error '%s'.\n", strerror(errno));
    _exit(-1);
  }
  waitpid(-1, (int*)&err, 0);
  err = WEXITSTATUS(err);

  if (err == 127)
    fprintf(stderr, "Failed to execute '%s'.\n", cmd), exit(1);

  if (fileExists(requested) == false)
    fprintf(stderr, "Failed to find or fetch file '%s'.\n", requested), exit(1);

  delete [] path;
  delete [] cmd;
}
