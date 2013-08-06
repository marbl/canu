////////////////////////////////////////
//
//  Library of reads.
//
//  You can add boolean flags and small-value types to the 64-bit-wide
//  bit-vectors immediately below.
//
//  IMPORTANT!  The default value of these should be 0.
//
//  This can make for some twisted logic, but consider adding an
//  on/off flag like doNotQVTrim and goodBadQVThreshold do -- if
//  doNotQVTrim is set, the value of goodBadQVThreshold is ignored,
//  otherwise, it is used.  Further, goodBadQVThreshold == 0 is
//  interpreted as the 'default' value (set in the client code),
//
//  IMPORTANT!  Try not to distrub the order.  It should be OK to add
//  new fields by stealing bits from the 'spare' field, and creating
//  the new field just after the spare.  The existing named fields
//  should be left in the same order and location.  (Though, adding
//  fields and remaining compatible with an existing store has not
//  actually been tested.)
//
//  When adding new fields, please update:
//    AS_GKP/AS_GKP_edit.C      -- to add edit interfaces
//    AS_PER/AS_PER_gkLibrary.C -- to add encoding/decoding of features
//    AS_PER/gkLibrary.H        -- this file
//
//  You need to also update clients so they know about the new
//  feature.  In general, only clients that act on a feature need to
//  know about it.
//

#define LIBRARY_NAME_SIZE 128
#define LIBRARY_NAME_DEFAULT "UNDEFINED"

class gkLibrary {
public:
  gkLibrary() {
    libraryUID = AS_UID_undefined();

    memset(libraryName, 0, sizeof(char) * LIBRARY_NAME_SIZE);
    strcpy(libraryName, LIBRARY_NAME_DEFAULT);

    mean       = 0.0;
    stddev     = 0.0;

    gkLibrary_clearFeatures();
  };
  ~gkLibrary() {
  };

  //  AS_MSG reads protoIO, turns a library into a LibraryMesg, which
  //  has the features and values strings.
  //
  //  AS_GKP receives the LibraryMesg, checks sanity, and converts it
  //  into gkLibrary.  It uses the functions below to populate the
  //  GKLR.
  //
  void gkLibrary_decodeFeatures(LibraryMesg *lmesg);
  void gkLibrary_encodeFeatures(LibraryMesg *lmesg);
  void gkLibrary_encodeFeaturesCleanup(LibraryMesg *lmesg);

  void gkLibrary_clearFeatures(void) {
    spareUN3                   = 0;
    spareUN2                   = 0;
    spareUN1                   = 0;

    spareUTG                   = 0;
    forceBOGunitigger          = 0;
    isNotRandom                = 0;

    spareALN                   = 0;
    doNotTrustHomopolymerRuns  = 0;

    spareOBT                   = 0;
    doTrim_initialNone         = 0;
    doTrim_initialMerBased     = 0;
    doTrim_initialFlowBased    = 0;
    doTrim_initialQualityBased = 1;  //  NOTE, NOT ZERO!

    doRemoveDuplicateReads     = 0;

    doTrim_finalLargestCovered = 0;
    doTrim_finalEvidenceBased  = 1;  //  NOTE, NOT ZERO!
    doTrim_finalBestEdge       = 0;

    doRemoveSpurReads          = 1;  //  NOTE, NOT ZERO!
    doRemoveChimericReads      = 1;  //  NOTE, NOT ZERO!
    doCheckForSubReads         = 0;

    //chimeraFromLinker          = 0;
    //chimeraFromBiotin          = 0;

    doConsensusCorrection      = 0;

    spareGKP                   = 0;
    forceShortReadFormat       = 0;

    spareLIB                   = 0;
    readsAreReversed           = 0;  //  NOT A USER SETTABLE FEATURE
    constantInsertSize         = 0;
    orientation                = 0;
  };

  void gkLibrary_setFeature(char *fea, char *val);

private:
public:
  AS_UID         libraryUID;
  char           libraryName[LIBRARY_NAME_SIZE];

  double         mean;
  double         stddev;

  //  Spare bits for new modules, one word for each
  uint64         spareUN3;
  uint64         spareUN2;
  uint64         spareUN1;

  //  Unitigger options
  uint64         spareUTG:62;                  //  Spare for unitigger changes
  uint64         forceBOGunitigger:1;          //  Default 0 -> allow any unitigger
  uint64         isNotRandom:1;                //  Default 0 -> fragments in here are random

  //  Alignment options (aka, consensus, but also overlapper)
  uint64         spareALN:63;                  //  Spare for alignment changes
  uint64         doNotTrustHomopolymerRuns:1;  //  Default 0 -> trust 'em

  //  OBT options
  uint64         spareOBT:52;                  //  Spare for OBT changes

  uint64         doTrim_initialNone:1;         //  Don't trim, use whatever clear range is present
  uint64         doTrim_initialMerBased:1;     //  For Illumina reads
  uint64         doTrim_initialFlowBased:1;    //  For 454 reads (NOT IMPLEMENTED)
  uint64         doTrim_initialQualityBased:1; //  For Sanger reads

  uint64         doRemoveDuplicateReads:1;     //  

  uint64         doTrim_finalLargestCovered:1; //  A simple 'pick the largest covered by overlaps'
  uint64         doTrim_finalEvidenceBased:1;  //  Pick the most plausible end points based on overlaps
  uint64         doTrim_finalBestEdge:1;

  uint64         doRemoveSpurReads:1;          //  
  uint64         doRemoveChimericReads:1;      //  
  uint64         doCheckForSubReads:1;         //  Default 0 -> reads don't contain subreads.  Added in v10

  uint64         doConsensusCorrection:1;      //  Default 0 -> fragments should be used with default correction

  //  Gatekeeper options
  uint64         spareGKP:63;                  //  Spare for gatekeeper changes
  uint64         forceShortReadFormat:1;       //  Default 0 -> store as long reads

  //  Library options
  uint64         spareLIB:59;                  //  Spare for library changes
  uint64         readsAreReversed:1;           //  Default 0 -> gatekeeper dump fastq as is
  uint64         constantInsertSize:1;         //  Default 0 -> allow changes to insert size
  uint64         orientation:3;                //  Default 0 -> AS_READ_ORIENT_UNKNOWN
};

