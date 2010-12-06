#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "sim4polishFile.H"


//  Global pointer used during construction of the index.
//  (polishRecordSortArray)
//
sim4polishFile::polishRecord *__prsa;

int
__prsaEST(const void *a, const void *b) {
  u32bit aid = __prsa[ *((u32bit*)a) ]._ESTiid;
  u32bit bid = __prsa[ *((u32bit*)b) ]._ESTiid;

  if (aid < bid) return(-1);
  if (aid > bid) return(1);
  return(0);
}

int
__prsaGEN(const void *a, const void *b) {
  u32bit aid = __prsa[ *((u32bit*)a) ]._GENiid;
  u32bit bid = __prsa[ *((u32bit*)b) ]._GENiid;

  if (aid < bid) return(-1);
  if (aid > bid) return(1);
  return(0);
}




sim4polishFile::sim4polishFile(char *path, sim4polishStyle style) {

  _path = new char [strlen(path) + 1];
  strcpy(_path, path);

  _file = new readBuffer(path);

  _style = style;

  _polishRecordLen = 0;
  _polishRecordMax = 0;
  _polishRecord    = 0L;
  _polishRecordEST = 0L;
  _polishRecordGEN = 0L;

  _maxEST = 0;
  _maxGEN = 0;
  _ESTiidLocation = 0L;
  _GENiidLocation = 0L;
}


sim4polishFile::~sim4polishFile() {
  delete [] _path;
  delete [] _polishRecord;
  delete [] _polishRecordEST;
  delete [] _polishRecordGEN;
  delete [] _ESTiidLocation;
  delete [] _GENiidLocation;
}


sim4polishList*
sim4polishFile::getEST(u32bit iid) {
  sim4polishList *l = new sim4polishList();

  if (iid >= _maxEST)
    //fprintf(stderr, "Invalid EST iid "u32bitFMT", max is "u32bitFMT"\n", iid, _maxEST), exit(1);
    return(l);

  sim4polish     *p = 0L;
  u32bit          i = _ESTiidLocation[iid];

  if (i != ~u32bitZERO) {
    setPosition(_polishRecordEST[i]);

    p = new sim4polish(_file, _style);

    while ((p) && (p->_numExons > 0) && (p->_estID == iid)) {
      l->push(p);
      i++;
      setPosition(_polishRecordEST[i]);
      p = new sim4polish(_file, _style);
    }
  
    delete p;
  }

  return(l);
}


sim4polishList*
sim4polishFile::getGEN(u32bit iid, u32bit lo, u32bit hi) {
  fprintf(stderr, "sim4polishFile::getGEN() not implemented.  Sorry.\n");
  exit(1);
  return(0L);
}


sim4polish*
sim4polishFile::getNext(void) {
  return(new sim4polish(_file, _style));
}


void
sim4polishFile::setPosition(u32bit ordinal) {

  if (_polishRecord == 0L)
    buildIndex();

  if (ordinal >= _polishRecordLen)
    fprintf(stderr, "Failed to reposition %s to record "u32bitFMT", only "u32bitFMT" records\n", _path, ordinal, _polishRecordLen), exit(1);

  _file->seek(_polishRecord[ordinal]._fileposition);
}


void
sim4polishFile::loadIndex(void) {
  char   magic[8] = {0};
  char   cigam[8] = { 's', '4', 'p', 'F', 'i', 'l', 'e', '1'};
  int    len = strlen(_path) + 32;
  char  *nam = new char [len];

  sprintf(nam, "%s.sim4polishFile", _path);

  if (fileExists(nam)) {
    errno = 0;
    FILE *F = fopen(nam, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", nam, strerror(errno)), exit(1);

    fread(&magic, sizeof(char), 8, F);
    if (strncmp(magic, cigam, 8) != 0)
      fprintf(stderr, "Failed to open '%s': Not a sim4polishFile!\n", nam), exit(1);

    fread(&_polishRecordLen, sizeof(u32bit), 1, F);

    _polishRecord    = new polishRecord [_polishRecordLen];
    _polishRecordEST = new u32bit       [_polishRecordLen];
    _polishRecordGEN = new u32bit       [_polishRecordLen];

    fread( _polishRecord,    sizeof(polishRecord), _polishRecordLen, F);
    fread( _polishRecordEST, sizeof(u32bit),       _polishRecordLen, F);
    fread( _polishRecordGEN, sizeof(u32bit),       _polishRecordLen, F);

    fread(&_maxEST, sizeof(u32bit), 1, F);
    fread(&_maxGEN, sizeof(u32bit), 1, F);

    _ESTiidLocation = new u32bit [_maxEST];
    _GENiidLocation = new u32bit [_maxGEN];

    fread( _ESTiidLocation, sizeof(u32bit), _maxEST, F);
    fread( _GENiidLocation, sizeof(u32bit), _maxGEN, F);

    if (errno)
      fprintf(stderr, "Failed to read '%s': %s\n", nam, strerror(errno)), exit(1);

    fclose(F);
  }

  delete [] nam;
}


void
sim4polishFile::saveIndex(void) {
  char   cigam[8] = { 's', '4', 'p', 'F', 'i', 'l', 'e', '1'};
  int    len = strlen(_path) + 32;
  char  *nam = new char [len];

  sprintf(nam, "%s.sim4polishFile", _path);

  errno = 0;
  FILE *F = fopen(nam, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", nam, strerror(errno)), exit(1);

  fwrite(&cigam, sizeof(char), 8, F);

  fwrite(&_polishRecordLen, sizeof(u32bit), 1, F);
  fwrite( _polishRecord,    sizeof(polishRecord), _polishRecordLen, F);
  fwrite( _polishRecordEST, sizeof(u32bit),       _polishRecordLen, F);
  fwrite( _polishRecordGEN, sizeof(u32bit),       _polishRecordLen, F);

  fwrite(&_maxEST, sizeof(u32bit), 1, F);
  fwrite(&_maxGEN, sizeof(u32bit), 1, F);
  fwrite( _ESTiidLocation, sizeof(u32bit), _maxEST, F);
  fwrite( _GENiidLocation, sizeof(u32bit), _maxGEN, F);

  if (errno)
    fprintf(stderr, "Failed to write '%s': %s\n", nam, strerror(errno)), exit(1);

  fclose(F);

  delete [] nam;
}


void
sim4polishFile::buildIndex(void) {

  loadIndex();

  if (_polishRecord == 0L) {
    fprintf(stderr, "sim4polishFile::buildIndex()-- building index for '%s'\n", _path);

    _file->seek(0);

    //  Allocate a bunch of space for stuff
    //
    _polishRecordLen = 0;
    _polishRecordMax = 3355443;  //  ~128MB for all three
    _polishRecord    = new polishRecord [_polishRecordMax];


    //  Read all polishes, storing stuff, reallocating more space if
    //  needed.
    //
    off_t       fp = _file->tell();
    sim4polish *p  = new sim4polish(_file, _style);

    while (p) {
      if (_polishRecordLen >= _polishRecordMax) {
        _polishRecordMax *= 2;
        polishRecord  *n = new polishRecord [_polishRecordMax];
        memcpy(n, _polishRecord, sizeof(polishRecord) * _polishRecordLen);
        delete [] _polishRecord;
        _polishRecord = n;
      }

      _polishRecord[_polishRecordLen]._fileposition = fp;
      _polishRecord[_polishRecordLen]._ESTiid = p->_estID;
      _polishRecord[_polishRecordLen]._GENiid = p->_genID;
      _polishRecord[_polishRecordLen]._GENlo  = p->_exons[0]._genFrom;
      _polishRecord[_polishRecordLen]._GENhi  = p->_exons[p->_numExons-1]._genTo;
      _polishRecordLen++;

      if ((_polishRecordLen & 0xfff) == 0) {
        fprintf(stderr, "polishes: "u32bitFMT"\r", _polishRecordLen);
        fflush(stderr);
      }

      delete p;
      
      fp = _file->tell();
      if (_file->eof())
        p = NULL;
      else 
        p  = new sim4polish(_file, _style);
    }


    //  Sort the indices by EST and GEN iid's.  Pain in the butt, we
    //  need to access _polishRecord to sort *EST and *GEN, but
    //  qsort() doesn't support that.
    //
    //  Three solutions:
    //  1) use a custom sort
    //  2) use a global pointer to _polishRecord
    //  3) use a temporary array holding the sort key and position
    //
    _polishRecordEST = new u32bit [_polishRecordLen];
    _polishRecordGEN = new u32bit [_polishRecordLen];

    for (u32bit i=0; i<_polishRecordLen; i++)
      _polishRecordEST[i] = _polishRecordGEN[i] = i;

    __prsa = _polishRecord;
    qsort(_polishRecordEST, _polishRecordLen, sizeof(u32bit), __prsaEST);
    qsort(_polishRecordGEN, _polishRecordLen, sizeof(u32bit), __prsaGEN);
    __prsa = 0L;


    //  Scan the sorted lists, record the first location of each iid
    //
    _maxEST = _polishRecord[ _polishRecordEST[_polishRecordLen-1] ]._ESTiid + 1;
    _maxGEN = _polishRecord[ _polishRecordGEN[_polishRecordLen-1] ]._GENiid + 1;
    _ESTiidLocation = new u32bit [_maxEST];
    _GENiidLocation = new u32bit [_maxGEN];

    for (u32bit i=0; i<_maxEST; i++)
      _ESTiidLocation[i] = ~u32bitZERO;
    for (u32bit i=0; i<_polishRecordLen; i++) {
      u32bit iid = _polishRecord[ _polishRecordEST[i] ]._ESTiid;
      if (_ESTiidLocation[iid] == ~u32bitZERO)
        _ESTiidLocation[iid] = i;
    }

    for (u32bit i=0; i<_maxGEN; i++)
      _GENiidLocation[i] = ~u32bitZERO;
    for (u32bit i=0; i<_polishRecordLen; i++) {
      u32bit iid = _polishRecord[ _polishRecordGEN[i] ]._GENiid;
      if (_GENiidLocation[iid] == ~u32bitZERO)
        _GENiidLocation[iid] = i;
    }


    //  Save the index
    //
    saveIndex();

    //  Be nice, reposition the file to the start.
    //
    _file->seek(0);
  }
}


