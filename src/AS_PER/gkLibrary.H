////////////////////////////////////////
//
//  Library of reads.
//
//  You can add boolean flags and small-value types to the
//  64-bit-wide bit-vector immediately below.  And just in case you
//  need A LOT of space, you've got ONE-HUNDRED-AND-TWENTY-EIGHT
//  bits!!
//
//  IMPORTANT!  The default value of these should be 0.  This can
//  make for some twisted logic, but consider adding an on/off flag
//  like doNotQVTrim and goodBadQVThreshold do.  The idea is that
//  value == 0 is interpreted as the 'default' value, and the on/off
//  restores what value 0 commonly means.
//
//  When adding new values, please update AS_PER_gkpStore.c,
//  AS_GKP_checkLibrary.c and sffToCA.c, thanks.
//
class gkLibrary {
public:
  gkLibrary() {};
  ~gkLibrary() {};

  //  AS_MSG reads protoIO, turns a library into a LibraryMesg, which
  //  has the features and values strings.
  //
  //  AS_GKP receives the LibraryMesg, checks sanity, and converts it
  //  into gkLibrary.  It uses the functions below to populate the GKLR.
  //
  void gkLibrary_decodeFeatures(LibraryMesg *lmesg);
  void gkLibrary_encodeFeatures(LibraryMesg *lmesg);
  void gkLibrary_encodeFeaturesCleanup(LibraryMesg *lmesg);

private:
public:
  AS_UID         libraryUID;

  double         mean;
  double         stddev;

  uint64         spare2:64;
  uint64         spare1:49;

  //  Unitigger options
  uint64         forceBOGunitigger:1;          //  Default 0 -> allow any unitigger
  uint64         isNotRandom:1;                //  Default 0 -> fragments in here are random

  //  Consensus and alignment options
  uint64         doNotTrustHomopolymerRuns:1;  //  Default 0 -> trust 'em

  //  OBT options
  uint64         doNotQVTrim:1;                //  Default 0 -> do initial qv trimming in OBT
  uint64         goodBadQVThreshold:6;         //  Default 0 -> threshold is 12 in initialTrim.C
  uint64         doNotOverlapTrim:1;           //  Default 0 -> do trimming

  //  Gatekeeper options
  uint64         useShortFragments:1;          //  Default 0 -> don't store gkShortFragment types

  //  Mate pairs
  uint64         orientation:3;                //  Default 0 -> AS_READ_ORIENT_UNKNOWN
};

