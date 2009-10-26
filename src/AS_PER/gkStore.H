
////////////////////////////////////////
//
//  The gkpStore
//
//  gkStoreInfo is saved on disk.
//  gkStore is the in memory structure used to access the data.
//
class gkStoreInfo {
private:
public:
  uint64    gkMagic;
  uint64    gkVersion;

  uint32    gkLibrarySize;
  uint32    gkPackedFragmentSize;
  uint32    gkNormalFragmentSize;
  uint32    gkStrobeFragmentSize;
  uint32    gkPlacementSize;

  uint32    readMaxLenBits;

  //  Statistics on our load

  uint32    libInput;
  uint32    libLoaded;
  uint32    libErrors;
  uint32    libWarnings;

  uint32    frgInput;
  uint32    frgLoaded;
  uint32    frgErrors;
  uint32    frgWarnings;

  uint32    lkgInput;
  uint32    lkgLoaded;
  uint32    lkgErrors;
  uint32    lkgWarnings;
   
  uint32    sffInput;
  uint32    sffLoaded;
  uint32    sffErrors;
  uint32    sffWarnings;

  uint32    sffLibCreated;

  uint32    plcInput;
  uint32    plcLoaded;
  uint32    plcErrors;
  uint32    plcWarnings;

  //  Counts of types of things we have loaded

  uint32    numRandom;

  //  The IID space is broken into three classes.  IIDs do NOT track
  //  the order fragments are loaded.  All short reads are first, then
  //  medium, then long.

  uint32    numPacked;
  uint32    numNormal;
  uint32    numStrobe;

  friend class gkFragment;
  friend class gkStore;
  friend class gkStream;
};


class gkStore {
private:
  gkStore(const char *path, int partnum);

  void         gkStore_open(int writable);
  void         gkStore_create(void);

public:
  gkStore();
  gkStore(const char *path, int creatable, int writable);
  ~gkStore();

  void         gkStore_clear(void);

  void         gkStore_load(AS_IID firstElem, AS_IID lastElem, int flags);
  void         gkStore_loadPartition(uint32 partition);

  void         gkStore_delete(void);

  const char  *gkStore_path(void) { return(storePath); };

  ////////////////////////////////////////
  //
  //
  //
  ////////////////////////////////////////

  void         gkStore_buildPartitions(short *partitionMap, uint32 maxPart);

  ////////////////////////////////////////
  //
  //  AS_PER_gkStore.c
  //
  ////////////////////////////////////////

  int32      gkStore_getNumLibraries(void)       { return(inf.libLoaded); }
  int32      gkStore_getNumFragments(void)       { return(inf.frgLoaded); }
  int32      gkStore_getNumRandomFragments(void) { return(inf.numRandom); }

  gkLibrary *gkStore_getLibrary(int id)                { return((id == 0) ? NULL : (gkLibrary *)getIndexStorePtr(lib, id)); }
  void       gkStore_getLibrary(int id, gkLibrary *lb) { getIndexStore(lib, id, lb); };
  void       gkStore_setLibrary(int id, gkLibrary *lb) { setIndexStore(lib, id, lb); }
  void       gkStore_addLibrary(AS_UID u, gkLibrary *lb) {
    appendIndexStore(lib, lb);
    gkStore_setUIDtoIID(u, getLastElemStore(lib), AS_IID_LIB);
    inf.libLoaded++;
  };

  void      gkStore_decodeTypeFromIID(AS_IID iid, uint32& type, uint32& tiid);

  void      gkStore_getFragmentData(gkStream *gst, gkFragment *fr, uint32 flags);

  void      gkStore_getFragment(AS_IID iid, gkFragment *fr, int32 flags);
  void      gkStore_setFragment(gkFragment *fr);
  void      gkStore_delFragment(AS_IID iid, bool deleteMateFrag=false);
  void      gkStore_addFragment(gkFragment *fr);

  ////////////////////////////////////////
  //
  //  Placement constraint support.
  //
  ////////////////////////////////////////
private:
   void        gkStore_loadFRGtoPLC(void);
   int         gkStore_setFRGtoPLC(AS_IID iid, int plid);

public:
  int32        gkStore_getNumPlacement(void)                 { return(inf.plcLoaded); }
  int          gkStore_getFRGtoPLC(AS_IID iid);
  
  gkPlacement *gkStore_getPlacement(int id)                  { return((id == 0) ? NULL : (gkPlacement *)getIndexStorePtr(plc, id)); }
  void         gkStore_getPlacement(int id, gkPlacement *pl) { getIndexStore(plc, id, pl); }  
  gkPlacement *gkStore_getReadPlacement(AS_IID iid)          { 
    int index = gkStore_getFRGtoPLC(iid);
    return (index == 0 ? NULL : gkStore_getPlacement(index));
  }
  void         gkStore_setPlacement(int id, gkPlacement *pl) { setIndexStore(plc, id, pl); }
  void         gkStore_addPlacement(gkPlacement *pl);

  ////////////////////////////////////////
  //
  //  Clear range support.
  //
  ////////////////////////////////////////

  //  Add a 'new' clear range.  Initial values are set to the version
  //  in the fragment.  Without enabling a range, only
  //  getClearRegion() works, and returns [1,0].
  //
  void      gkStore_enableClearRange(uint32 which);
  void      gkStore_purgeClearRange(uint32 which);

  ////////////////////////////////////////
  //
  //  AS_PER_gkStore_UID.c
  //
  ////////////////////////////////////////

  //  Privates to load/save data
private:
  void           gkStore_loadUIDtoIID(void);
  void           gkStore_loadIIDtoUID(void);
  void           gkStore_loadSTRtoUID(void);


  //  These for general consumption
public:
  AS_IID         gkStore_getUIDtoIID(AS_UID uid, uint32 *type);
  int            gkStore_setUIDtoIID(AS_UID uid, AS_IID iid, uint32 type);
  AS_UID         gkStore_getIIDtoUID(AS_IID iid, uint32 type);

  //  This for gatekeeper the program
  void           gkStore_rebuildUIDtoIID(void);

  //  The next three functions should only be used by AS_UTL_UID.c (and the null constructor).
  AS_UID         gkStore_getUIDfromString(char *uidstr);
  char          *gkStore_getUIDstring(AS_UID uid);
  AS_UID         gkStore_addUID(char *uidstr);

  ////////////////////////////////////////
  //
  //
  //
  ////////////////////////////////////////

  void   gkStore_flagWarning(char t);
  void   gkStore_flagError(char t);

private:
  char                     storePath[FILENAME_MAX];

  uint32                   isReadOnly;
  uint32                   isCreating;

public:
  //  Sigh, needed for AS_GKP.
  gkStoreInfo              inf;

private:
  StoreStruct             *fpk;  //  Small fragments

  StoreStruct             *fnm;  //  Medium size fragments
  StoreStruct             *snm;
  StoreStruct             *qnm;

  StoreStruct             *fsb;  //  Large fragments
  StoreStruct             *ssb;
  StoreStruct             *qsb;

  StoreStruct             *lib;

  StoreStruct             *uid;
  
  // Stores placement constraints for fragments
  StoreStruct             *plc;
  HashTable_AS            *FRGtoPLC;

  //  Maps UIDs to IIDs; maps strings to UIDs.  The STR array
  //  holds ALL the string UIDs.
  //
  HashTable_AS            *UIDtoIID;
  HashTable_AS            *STRtoUID;

  //  We can generate a quick mapping from IID to UID for fragments.
  //
  uint64                  *frgUID;

  //  Clear ranges
  //
  gkClearRange           **clearRange;

  //  These are for use when building a gkpStore.  They map an IID
  //  to a type and a tiid.
  //
  uint32                   IIDmax;
  uint8                   *IIDtoTYPE;
  uint32                  *IIDtoTIID;

  //  The rest are for a partitioned fragment store.
  //
  //  We load all frg and qlt in this partition into memory.  The map
  //  converts an iid (global fragment iid) into a pointer to the
  //  correct frg record, which we can then use to grab the
  //  encoded seq/qlt.
  //
  int32                    partnum;

  StoreStruct             *partfpk;
  StoreStruct             *partfnm, *partqnm;
  StoreStruct             *partfsb, *partqsb;

  HashTable_AS            *partmap;

  friend class gkStream;
  friend class gkFragment;  //  for clearRange
};
