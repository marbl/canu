
seqFactory *seqFactory::me = 0L;


seqFactory::seqFactory() {
  _filesNum = 0;
  _filesMax = 16;
  _files = new seqFile * [_filesMax];

  registerFile(new fastaFile);
}


seqFactory::~seqFactory() {
  for (u32bit i=0; i<_filesNum; i++)
    delete _files[i];
  delete [] _files;
}


void           
seqFactory::registerFile(seqFile *f) {
  if (_filesNum >= _filesMax) {
    fprintf(stderr, "Hmmmm!  Wow!  You registered lots of files!  Now fix %s at line %d.\n", __FILE__, __LINE__);
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

  fprintf(stderr, "ERROR: Cannot determine type of file '%s'.\n", name);
  exit(1);
  return(n);
}
