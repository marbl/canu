#ifndef SIM4POLISHFILE
#define SIM4POLISHFILE

#include "sim4polish.H"
#include "sim4polishList.H"

//
//  (original motivation)
//  Needed: something to return polishes from a file without reading
//  in everything, and without doing lots of I/O.
//
//    polishFile   A(name)
//    sim4polishList  l = A.getEST(id);
//    sim4polishList  l = A.getGEN(id, lo=0, hi=0);
//
//    sim4polish      x = A.getNext();
//
//    A.seek(polish-ordinal to seek to);
//
//    On first invocation, it reads the whole file, building a map of
//    ESTid, GENid to file position.  Maybe also store GENlo and GENhi
//    in this map.  Map is stored on disk as "file.polishFileMap"
//


//
//  The first call to getEST(), getGEN() or seek() will cause a
//  while file scan to be performed.  From this, we build a list of
//  all the polishes, and their iid's.  This index is cached on disk
//  as '${path}.sim4polishFile'
//

class sim4polishFile {
public:
  sim4polishFile(char *path, sim4polishStyle style);
  ~sim4polishFile();

  //  Number of EST's actually present
  //
  uint32    maxIID(void) {
    return(_maxEST);
  };

  //  Returns all matches with:
  //    ESTid == iid
  //    GENid == iid AND that overlap the range lo...hi
  //
  //  N.B. getNext() doesn't really mean much after these.
  //
  sim4polishList  *getEST(uint32 iid);
  sim4polishList  *getGEN(uint32 iid, uint32 lo=0, uint32 hi=~uint32ZERO);

  //  Returns the next polish in the file.
  //
  sim4polish      *getNext(void);

  //  Positions the file pointer to the ordinal'th polish in the file.
  //  getNext() will then return that polish.
  //
  void             setPosition(uint32 ordinal);

  //  24 bytes/record.  A typical large EST set has ~10 million
  //  matches, this fits into 240MB.  That's reasonable.
  //
  //  I really hate to make this public, but sorting needs it.
  //
  struct polishRecord {
    off_t      _fileposition;
    uint32     _ESTiid;
    uint32     _GENiid;
    uint32     _GENlo;
    uint32     _GENhi;
  };

private:
  char            *_path;
  readBuffer      *_file;

  //  One record for each polish, in the order they are in the file.
  //  One integer pointer into _polishRecord, sorted by either the EST
  //  or GEN iid.
  //
  uint32           _polishRecordLen;
  uint32           _polishRecordMax;
  polishRecord    *_polishRecord;
  uint32          *_polishRecordEST;
  uint32          *_polishRecordGEN;
  
  sim4polishStyle  _style;

  //  One integer pointer for each iid we've seen.  Pointer into
  //  _polishRecordEST or _polishRecordGEN.  If memory is tight, we
  //  could binary search those arrays instead.
  //
  uint32           _maxEST;
  uint32           _maxGEN;
  uint32          *_ESTiidLocation;
  uint32          *_GENiidLocation;

  void             loadIndex(void);
  void             saveIndex(void);
  void             buildIndex(void);
};



#endif  //  SIM4POLISHFILE
