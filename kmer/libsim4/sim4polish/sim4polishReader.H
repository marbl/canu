#ifndef SIM4POLISHREADER
#define SIM4POLISHREADER

#include "sim4polish.H"
#include "util++.H"

//  Simple class to read the contents of a file of alignments.  The file can be either sim4db, gff3
//  or atac format.  No support for random access is provided, just sequential access.

class sim4polishWriter;

class sim4polishReader {
public:
  sim4polishReader(const char *name, sim4polishWriter *writer=0L);
  ~sim4polishReader();

  //  Returns the next alignment in the file.  NULL is returned if there are no more alignments.
  //
  sim4polish *nextAlignment(void);
  bool        nextAlignment(sim4polish * &p);

  sim4polishStyle getsim4polishStyle(void) { return _style; } 

private:
  readBuffer      *_rb;
  sim4polishStyle  _style;
};

#endif // SIM4POLISHREADER
