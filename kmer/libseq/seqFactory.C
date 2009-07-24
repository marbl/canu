#include "seqFactory.H"

#include "fastaStdin.H"
#include "fastaFile.H"
#include "seqStore.H"

seqFactory *seqFactory::me = 0L;


seqFactory::seqFactory() {
  _filesNum = 0;
  _filesMax = 16;
  _files = new seqFile * [_filesMax];

  registerFile(new fastaStdin);
  registerFile(new fastaFile);
  registerFile(new seqStore);
  //registerFile(new sffFile);
}


seqFactory::~seqFactory() {
  for (u32bit i=0; i<_filesNum; i++)
    delete _files[i];
  delete [] _files;
}


void           
seqFactory::registerFile(seqFile *f) {
  if (_filesNum >= _filesMax) {
    fprintf(stderr, "seqFactory::registerFile()--  Wow!  You registered lots of files!  Now fix %s at line %d.\n", __FILE__, __LINE__);
    exit(1);
  }
  _files[_filesNum++] = f;
}


seqFile *
seqFactory::openFile(const char *name) {
  seqFile  *n = 0L;

  for (u32bit i=0; i<_filesNum; i++) {
    n = _files[i]->openFile(name);
    if (n)
      return(n);
  }

  fprintf(stderr, "seqFactory::registerFile()--  Cannot determine type of file '%s'.  Tried:\n", name);

  for (u32bit i=0; i<_filesNum; i++)
    fprintf(stderr, "seqFactory::registerFile()--         '%s'\n", _files[i]->getFileTypeName());

  exit(1);
  return(n);
}
