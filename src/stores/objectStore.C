
/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "types.H"
#include "arrays.H"
#include "strings.H"

#include "objectStore.H"

#include <libgen.h>
#include <sys/wait.h>



extern char **environ;  //  Where, oh where, is this really defined?!



static
char *
findSeqStorePath(char *requested) {
  splitToWords  F(requested, splitPaths);

  if (F.numWords() < 2)
    return(NULL);

  char  *filename  = F.last(0);
  char  *storename = F.last(1);

  //  If not a blobs file name, return no file.

  if ((strlen(filename) != 10) ||
      (filename[0] != 'b') ||
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

  //  If not an overlap store data file name, return no file.
  //  This format is defined in ovFile::createDataName() in ovStoreFile.C

  if ((strlen(filename) != 8) ||
      (isdigit(filename[0]) == 0) ||
      (isdigit(filename[1]) == 0) ||
      (isdigit(filename[2]) == 0) ||
      (isdigit(filename[3]) == 0) ||
      (filename[4] != '-')        ||
      (isdigit(filename[5]) == 0) ||
      (isdigit(filename[6]) == 0) ||
      (isdigit(filename[7]) == 0))
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



bool
fetchFromObjectStore(char *requested) {

  //  Decide if we even need to bother.  If the file exists locally, or if
  //  one of the environment variables is missing, no, we don't need to bother.

  if (fileExists(requested))
    return(false);

  char  *da = getenv("CANU_OBJECT_STORE_CLIENT_DA");
  char  *ns = getenv("CANU_OBJECT_STORE_NAMESPACE");
  char  *pr = getenv("CANU_OBJECT_STORE_PROJECT");

  if ((da == NULL) ||
      (ns == NULL) ||
      (pr == NULL))
    return(false);

  //  Try to figure out the object store path for this object based on the name
  //  of the requested file.  Paths to stores are relative, but we need them
  //  rooted in the assembly root directory:
  //
  //      ../../asm.seqStore -> ./asm.seqStore
  //      ../asm.ovlStore    -> ./correction/asm.ovlStore
  //
  //  For the seqStore, we can just grab the last two components.
  //  For the ovlStore, we need to parse out the subdirectory the store is in.

  char *path   = NULL;

  if (path == NULL)
    path = findSeqStorePath(requested);

  if (path == NULL)
    path = findOvlStorePath(requested);

  if (path == NULL)
    fprintf(stderr, "fetchFromObjectStore()-- requested file '%s', but don't know where that is.\n", requested), exit(1);

  //  With the path to the object figured out, finish making the path by appending
  //  the PROJEXT and NAMESPACE.

  char *object = new char [FILENAME_MAX+1];

  snprintf(object, FILENAME_MAX, "%s:%s/%s", pr, ns, path);

  //  Then report what's going on.

  fprintf(stderr, "fetchFromObjectStore()-- fetching file '%s'\n", requested);
  fprintf(stderr, "fetchFromObjectStore()--   from object '%s'\n", object);

  //  Build up a command we can execute after forking.

  char *args[8];

  args[0] = basename(da);
  args[1] = duplicateString("download");        //  Thanks, execve, for wanting mutable
  args[2] = duplicateString("--overwrite");     //  strings and making us jump through
  args[3] = duplicateString("--no-progress");   //  a hoop to get them without compiler
  args[4] = duplicateString("--output");        //  warnings.
  args[5] = requested;
  args[6] = object;
  args[7] = NULL;

  //  Fork and run the child command if we're the child.  Normally, evecve()
  //  doesn't return (because it obliterated the process it could return to).
  //  If it does return, an error occurred, so we just go BOOM too.  As per
  //  the manpage, _exit() MUST be used instead of exit(), so that
  //  stdin/out/err are left intact.
  //
  //  vfork() is dangerous.  If we're the child, all we're allowed to do
  //  after the call is execve() or _exit().  Absolutely nothing else.

  pid_t pid = vfork();

  if (pid == 0) {
    execve(da, args, environ);
    fprintf(stderr, "fetchFromObjectStore()-- execve() failed with error '%s'.\n", strerror(errno));
    _exit(127);
  }

  if (pid == -1)
    fprintf(stderr, "fetchFromObjectStore()-- vfork() failed with error '%s'.\n", strerror(errno)), exit(1);

  //  Otherwise, we're still the parent; wait for the child process to
  //  terminate.

  int   status = 0;
  pid_t wid    = waitpid(pid, &status, 0);

  if (wid == -1)
    fprintf(stderr, "fetchFromObjectStore()-- waitpid() failed with error '%s'.\n", strerror(errno)), exit(1);

  if ((WIFEXITED(status)) &&
      (WEXITSTATUS(status) == 127))
    fprintf(stderr, "fetchFromObjectStore()-- execve() failed to run the command.\n"), exit(1);

  //  If no file, it's fatal.
  if (fileExists(requested) == false)
    fprintf(stderr, "fetchFromObjectStore()-- failed fetch file '%s'.\n", requested), exit(1);

  delete [] args[1];
  delete [] args[2];
  delete [] args[3];
  delete [] args[4];

  delete [] path;
  delete [] object;

  return(true);
}
