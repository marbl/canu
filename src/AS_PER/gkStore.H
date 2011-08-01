
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
  uint32    gkPackedSequenceSize;
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

  //  The IID space is broken into three classes.  See comments in AS_PER_gkStore_IID.C

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

  void         gkStore_open(int writable) { gkStore_open(writable, 1); }
  void         gkStore_open(int writable, int doNotUseUIDs);
  void         gkStore_create(void);

  void         gkStore_clear(void);
  
  void         gkStore_construct(const char * path, int creatable, int writable, int doNotUseUIDs);

public:
  gkStore();
  gkStore(const char *path, int creatable, int writable);
  gkStore(const char *path, int creatable, int writable, int doNotUseUIDs);
  ~gkStore();

  //  The 'name' of this gkpStore.
  const char  *gkStore_path(void) { return(storePath); };

  ////////////////////////////////////////
  //
  //  Operatons on the whole store.  Loading into memory, removing from disk.
  //
  ////////////////////////////////////////

  void         gkStore_load(AS_IID firstElem, AS_IID lastElem, int flags);

  void         gkStore_loadPartition(uint32 partition);
  void         gkStore_buildPartitions(short *partitionMap, uint32 maxPart);

  void         gkStore_delete(void);

  uint64       gkStore_metadataSize(void);
  void         gkStore_metadataCaching(bool enable=true);

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
  void      gkStore_addIIDtoTypeMap(AS_IID iid, uint32 type, uint32 tiid);
  void      gkStore_computeRanges(AS_IID  bgnIID, AS_IID  endIID,
                                  uint32 &bgnPK,  uint32 &endPK,  uint32 &valPK,
                                  uint32 &bgnNM,  uint32 &endNM,  uint32 &valNM,
                                  uint32 &bgnSB,  uint32 &endSB,  uint32 &valSB);

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

public:    //  Sigh, public needed for AS_GKP.
  gkStoreInfo              inf;

private:
  StoreStruct             *fpk;  //  Packed fragments
  StoreStruct             *qpk;

  StoreStruct             *fnm;  //  Normal fragments
  StoreStruct             *snm;
  StoreStruct             *qnm;

  StoreStruct             *fsb;  //  Strobe fragments
  StoreStruct             *ssb;
  StoreStruct             *qsb;

  StoreStruct             *lib;

  StoreStruct             *uid;
  
  // Stores placement constraints for fragments
  StoreStruct             *plc;
  HashTable_AS            *FRGtoPLC;

  //  Maps UIDs to IIDs; maps strings to UIDs.  The STR array holds ALL the string UIDs.
  //
  HashTable_AS            *UIDtoIID;
  HashTable_AS            *STRtoUID;

  //  We can generate a quick mapping from IID to UID for fragments.
  //
  uint64                  *frgUID;

  //  Flag to indicate whether we are using UIDs. If this is set to true, UIDs will not be loaded
  uint32                   doNotLoadUIDs;
  
  //  Clear ranges
  //
  gkClearRange           **clearRange;

  //  If this is defined (not NULL) then the store contains reads of mixed types, AND, it was
  //  created 'out of order' (in order: packed first, then normal, then strobe).  See comments in
  //  AS_PER_gkStore_IID.C
  //
  uint32                   IIDmax;
  uint8                   *IIDtoTYPE;
  uint32                  *IIDtoTIID;

  //  The rest are for a partitioned fragment store.
  //
  //  We load all frg and qlt in this partition into memory.  The map converts an iid (global
  //  fragment iid) into a pointer to the correct frg record, which we can then use to grab the
  //  encoded seq/qlt.
  //
  int32                    partnum;

  StoreStruct             *partfpk, *partqpk;
  StoreStruct             *partfnm, *partqnm;
  StoreStruct             *partfsb, *partqsb;

  HashTable_AS            *partmap;

  friend class gkStream;
  friend class gkFragment;  //  for clearRange
  friend class gkClearRange;
};
