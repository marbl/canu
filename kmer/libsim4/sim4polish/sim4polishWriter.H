#ifndef SIM4POLISHWRITER
#define SIM4POLISHWRITER

#include "sim4polish.H"

//  Simple class for writing a file of alignments.  The file can be either sim4db, gff3 or atac
//  format.  This class makes sure that the file has a header (if needed) and takes care of
//  generating unique IDs for each gff3 file.

class sim4polishWriter {
public:
  sim4polishWriter(const char *name, sim4polishStyle style, bool hidden=false);
  ~sim4polishWriter();

private:
  //  If this was opened as a 'hidden' file, this is the only way to ever see the results again.
  //  The destructor complains if this method is never called on a hidden file.
  //
  FILE  *surrenderToReader(void);

  friend class sim4polishReader;

public:
  //  The source name is listed in column 2 of a GFF3.  It is supposed to be the name of the
  //  program that generated these alignments.
  //
  //  The source name MUST be shorter than 32 letters, and if it contains /'s (a path) only the
  //  last component is used.
  //
  void   setSourceName(const char *sourceName);

  //  The match ID prefix is used in column 9, when constructing the file-unique ID for each
  //  alignment.  If not set, it will default to the sourceName + a short salt derived from the
  //  current time and process id + an integer count starting at zero.  The default is chosen so
  //  that the ID's in resulting files are more-or-less globally unique.
  //
  //  The match ID prefix MUST be shorter than 32 letters.
  //
  void   setMatchIDPrefix(const char *prefix);

  //  Add an alignment to the file.
  //
  void   writeAlignment(sim4polish *out);

private:
  char             _otName[FILENAME_MAX];
  FILE            *_otFile;

  sim4polishStyle  _style;

  char             _sourceName[32];

  char             _matchIDprefix[32];
  char             _matchIDsalt[8];
  u64bit           _matchID;
};

#endif // SIM4POLISHWRITER
