
////////////////////////////////////////
//
//  gkClearRanges
//
//  This manages a disk-based list of clear ranges.  All fragments are
//  in the list.
//
//  It is used directly by OBT and ECR.
//
//  The constructor defers all loads until they're used.
//  The descructor updates all disk-based stores.
//
//  Both get() and set() will load/create a clear range on demand.
//
//  1) If the clear range isn't loaded, this will automagically load
//  it for us.
//
//  2) If there isn't a clear range to load, one will be created and
//  initialized to the latest, if we're allowed to create.  If we're
//  not allowed to create, set() fails (asserts) here.  Enable creation
//  with gkStore_enableClearRange(name), which will call
//  gkClearRange_enableCreate().
//
//  3) If not allowed to create, then an invalid clear range is
//  returned.
//

class gkClearRange {
public:
  gkClearRange(gkStore *gkp_, uint32 clearType_, uint32 create_=FALSE);
  ~gkClearRange();

  void       gkClearRange_purge(void);

  void       gkClearRange_enableCreate(void) { create = TRUE; };

  void       gkClearRange_getClearRegion(gkFragment *fr, uint32& begin, uint32& end);
  void       gkClearRange_setClearRegion(gkFragment *fr, uint32  begin, uint32  end);

private:
  //  Load a range from disk, or, if create, copy from latest.
  void       gkClearRange_configurePacked(void);
  void       gkClearRange_configureNormal(void);
  void       gkClearRange_configureStrobe(void);

  //  These allocate space for new clear ranges (used in gatekeeper
  //  when loading fragments) and initialize them to [1,0].  The dirty
  //  flag is NOT set, so unless someone actually sets one of the
  //  ranges, these empty ranges are not saved.
  //
  void       gkClearRange_makeSpacePacked(AS_IID tiid, uint32 bgn, uint32 end);
  void       gkClearRange_makeSpaceNormal(AS_IID tiid, uint32 bgn, uint32 end);
  void       gkClearRange_makeSpaceStrobe(AS_IID tiid, uint32 bgn, uint32 end);

  gkStore   *gkp;

  uint32     clearType;
  uint32     create;

  //  configured -- 1 if we've tried to load, or have loaded.
  //  dirty      -- 1 if we need to update the disk
  //  maxiid     -- the largest iid we have space for
  //                we allocate 2*max + 2 
  //
  uint32     pkconfigured;
  uint32     pkdirty;
  uint64     pkmaxiid;
  uint8     *pk;

  uint32     nmconfigured;
  uint32     nmdirty;
  uint64     nmmaxiid;
  uint16    *nm;

  uint32     sbconfigured;
  uint32     sbdirty;
  uint64     sbmaxiid;
  uint32    *sb;

  friend class gkStore;
};

