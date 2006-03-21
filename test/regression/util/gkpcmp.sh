######################################################################
# Compare the gatekeeper store files created by two separate
# gatekeeper instances when they are fed the same input frag
# file. Comparison is done using the md5sum of the store files.
######################################################################
#!/bin/sh

###
# Name by which this script was invoked.
###
progname=`echo "$0" | sed -e 's/[^\/]*\///g'`

usage="Usage: $progname -f <frag> -a <gkpbin> -b <gkpbin>

Options are:
-f <file>            Specify the frag file to run against
-a <file>            Specify 1st gatekeeper binary
-b <file>            Specify 2nd gatekeeper binary
"

FRAG_PREF=""
FRAG_FILE=""

RUN1_MD5S=/tmp/run1_md5sums.$$
RUN2_MD5S=/tmp/run2_md5sums.$$

XXD_BIN=xxd

###
# $GKP_RUN1_BIN   1st gatekeeper binary
# $GKP_RUN2_BIN   2nd gatekeeper binary
###
GKP_RUN1_BIN=""
GKP_RUN2_BIN=""

CLEANUP() {
  rm -rf RUN1_${FRAG_PREF}.gkpStore
  rm -rf RUN2_${FRAG_PREF}.gkpStore

  rm -f  ${FRAG_PREF}.err
  rm -f  ${FRAG_PREF}.ign
  rm -f  ${FRAG_PREF}.inp

  rm -f  $RUN1_MD5S
  rm -f  $RUN2_MD5S

  file $FRAG_FILE | grep "symbolic link" > /dev/null 2>&1
  ret=$?
  if [ "$ret" = "0" ]; then
    rm -f $FRAG_FILE
  fi
  echo
}

PASSED() {
  echo "-----------------------------"; 
  echo "*** SUCCESS - Passed Test ***"; 
  echo "-----------------------------"; echo ""

  CLEANUP

  exit 0
}

FAILED() {
  echo "-----------------------------"; 
  echo "*** FAILURE - Failed Test ***"; 
  echo "-----------------------------"; echo ""

  ###
  # no cleanup so user can examine files
  # for cause of error
  ###
  exit 1
}

#####
# Function 
#   patch_gkpStore_tfields
# 
# Purpose 
#   The gatekeeper store files contain creation
#   time and last update time fields within their 
#   headers which makes files that are otherwise the 
#   same have different md5 sums. So, we go through 
#   and zero out the time fields of these files 
#   before running md5sum against them
# 
# Arguments
#   file - file to patch
#
# Return
#   success - NONE
#   failure - exit 1
#####
patch_gkpStore_tfields()
{
  local file=$1

  if [ ! -f $file ]; then
    echo "patch_gkpStore_tfields: $file does not exist"
    exit 1
  fi

  ####
  # The header for these store files follows. The
  # time field is 40 bytes into this header. This
  # is the field we zero out:
  #
  # 32-bit system layout
  # --------------------
  # typedef struct{                   # bytes   offset
  #  unsigned int isDeleted:1;        -------   ------
  #  unsigned int type:3;
  #  unsigned int :28;                // 4      0    0x00
  #  unsigned int :32;                // 4      4    0x04
  #  char storeType[8];               // 8      8    0x08
  #  int64 firstElem;                 // 8      16   0x10
  #  int64 lastElem;                  // 8      24   0x18
  #  int32 version;                   // 4      32   0x20
  #  int32 elementSize;               // 4      36   0x24
  #  time_t creationTime;             // 4      40   0x28
  # #ifndef __x86_64__                
  #  uint32 padtime_t1;               // 4      44   0x2c
  # #endif
  #  time_t lastUpdateTime;           // 4      48   0x30
  # #ifndef __x86_64__
  #  uint32 padtime_t2;               // 4      52   0x34
  #endif

  ####
  echo '0000028: 00000000' | $XXD_BIN -r - $file
  echo '0000030: 00000000' | $XXD_BIN -r - $file

  ####
  # Patch the padding fields: if we're on a 32-bit box
  # then this is redundant. If we're on 64-bit box though
  # then this is required as our time fields are twice
  # as large; IE,
  #
  # typedef struct{
  #  ...
  #  time_t creationTime;             // 8      40   0x28
  #  time_t lastUpdateTime;           // 8      48   0x30
  ####
  echo '000002c: 00000000' | $XXD_BIN -r - $file
  echo '0000034: 00000000' | $XXD_BIN -r - $file
}

#####
# Function 
#   patch_gkpStore_phash
# 
# Purpose 
#   The gatekeeper store phash file's header
#   contains some fields which are pointers set
#   from the return value of [mc]alloc(). As such
#   their values can change from run to run. So,
#   we zero these fields out so that they do not
#   skew the md5sum
# 
# Arguments
#   NONE
#
# Return
#   success - NONE
#   failure - exit 1
#####
patch_gkpStore_phash()
{
  local file=$1

  if [ ! -f $file ]; then
    echo "patch_gkpStore_phash: $file does not exist"
    exit 1
  fi

  ####
  # The following fields have to be zero'd out:
  #   PHashNode_AS * allocated    @ offset 0x68
  #   int32        * buckets      @ offset 0x70
  #   char         * fileName     @ offset 0x78
  #   FILE         * fp           @ offset 0x80
  #
  #
  # 32-bit system layout
  # --------------------
  #                                       # bytes  offset
  #                                       -------  ------
  # typedef struct{              
  #   int32 numBuckets;                   // 4     0    0x00
  #   int32 freeList;                     // 4     4    0x04
  # 
  #   int32 numNodes;                     // 4     8    0x08
  #   int32 numNodesAllocated;            // 4     12   0x0c
  #  
  #   CDS_UID_t lastKey;                  // 8     16   0x10
  #  
  #   int32 lastNodeAllocated;            // 4     24   0x18
  #   int32 collisions;                   // 4     28   0x1c
  # 
  #   uint32 hashmask;                    // 4     32   0x20
  #   int32 dummy1For8byteWordPadding;    // 4     36   0x24
  # 
  #   CDS_IID_t counts[1<<LOG_NUM_TYPES]; // 64    40   0x28
  # 
  #   PHashNode_AS *allocated;            // 4     104  0x68
  #   void *dummyPadPtr1;                 // 4     108  0x6c
  # 
  #   int32 *buckets;                     // 4     112  0x70
  #   void *dummyPadPtr2;                 // 4     116  0x74
  # 
  #   char *fileName;                     // 4     120  0x78
  #   void *dummyPadPtr3;                 // 4     124  0x7c
  # 
  #   FILE *fp;                           // 4     128  0x80
  #   void *dummyPadPtr4;                 // 4     132  0x84
  # 
  #   int32 isDirty;                      // 4     136  0x88
  #   int32 isReadWrite;                  // 4     140  0x8c
  # } PHashTable_AS;
  ####
  echo '0000068: 00000000' | $XXD_BIN -r - $file
  echo '0000070: 00000000' | $XXD_BIN -r - $file
  echo '0000078: 00000000' | $XXD_BIN -r - $file
  echo '0000080: 00000000' | $XXD_BIN -r - $file

  ####
  # Patch the padding fields: if we're on a 32-bit box
  # then this is redundant. If we're on 64-bit box though
  # then this is required as our pointer fields are twice
  # as large; IE,
  #
  # typedef struct{
  #  ...
  #   PHashNode_AS *allocated;            // 8     104  0x68
  #   int32 *buckets;                     // 8     112  0x70
  #   char *fileName;                     // 8     120  0x78
  #   FILE *fp;                           // 8     128  0x80
  ####
  echo '000006c: 00000000' | $XXD_BIN -r - $file
  echo '0000074: 00000000' | $XXD_BIN -r - $file
  echo '000007c: 00000000' | $XXD_BIN -r - $file
  echo '0000084: 00000000' | $XXD_BIN -r - $file
}

#####
# Function
#   gen_md5s
#
# Purpose
#   Generate MD5 sums of gatekeeper created store
#   files by running 'md5sum' against the files 
#   contained in the given directory
# 
# Arguments
#   dir - directory containing gatekeeper store
#         files (*.gkpStore/)
#   out - output file name
#
# Return
#   success - NONE
#   failure - exit 1
#####
gen_md5s()
{
  local dir=$1
  local out=$2

  if [ ! -d $dir ]; then
    echo "gen_md5s: directory $dir does not exist";
    exit 1;
  fi

  pushd $dir > /dev/null 2>&1
  for file in $(ls -1 | grep -v "\.phash$")
  do
    patch_gkpStore_tfields $file
    md5sum $file >> $out
  done
  for file in $(ls -1 | grep "\.phash$")
  do
    patch_gkpStore_phash $file
    md5sum $file >> $out
  done
  popd > /dev/null 2>&1
}

#####
# Function
#   gen_run1_md5s
#
# Purpose
#   Generate the md5s for the first run of
#   gatekeeper against the input
# 
# Arguments
#   
# Return
#   success - NONE
#   failure - exit 1
#####
gen_run1_md5s()
{
  local dir=RUN1_${FRAG_PREF}.gkpStore
  local out=$RUN1_MD5S

  gen_md5s $dir $out
}

#####
# Function
#   gen_run2_md5s
#
# Purpose
#   Generate the md5s for the second run of
#   gatekeeper against the input
# 
# Arguments
#   
# Return
#   success - NONE
#   failure - exit 1
#####
gen_run2_md5s()
{
  local dir=RUN2_${FRAG_PREF}.gkpStore
  local out=$RUN2_MD5S

  gen_md5s $dir $out
}

#####
# Function
#   run_gkpbin
#
# Purpose
#   Run gatekeeper binary against the user 
#   provided frag file
# 
# Arguments
#   gkp - gatekeeper binary to use
#   run - run number
#
# Return
#   success - NONE
#   failure - exit 1
#####
run_gkpbin()
{
  local GKP_BIN=$1
  local RUN_NUM=$2

  if [ ! -f $GKP_BIN ]; then
    echo "run_gkpbin: $GKP_BIN executable not found"
    FAILED
  fi

  if [ -z "$RUN_NUM" ]; then
    echo "run_gkpbin: RUN_NUM not set"
    FAILED
  fi

  if [ ! -f $FRAG_FILE ]; then
    echo "run_gkpbin: $FRAG_FILE not found"
    FAILED
  fi

  ###
  # Run gatekeeper against the frag file provided
  # by the user
  ###
  echo "COMMAND: $GKP_BIN -X -C -N -Q -P -f ${FRAG_PREF}.gkpStore $FRAG_FILE"

  $GKP_BIN -X -C -N -Q -P -f \
    ${FRAG_PREF}.gkpStore $FRAG_FILE > /dev/null 2>&1


  ret=$?
  if [ "$ret" -ne "0" ]; then
    echo "run_gkpbin: $GKP_BIN failed!"
    exit 1
  fi

  if [ ! -d ${FRAG_PREF}.gkpStore ]; then
    echo "run_bkpbin: Frag store NOT created"
    exit 1
  fi
  
  ###
  # Move the gatekeeper store directory into a directory
  # prefixed with the run number identifying this run
  ###
  mv ${FRAG_PREF}.gkpStore RUN${RUN_NUM}_${FRAG_PREF}.gkpStore
}

#####
# Function
#   RUN_TEST
#
# Purpose
#   Run the given gatekeeper binary against
#   the frag file within this script. Verify
#   that the results match the baseline results
#   by comparing MD5 checksums of the generated
#   gatekeeper store files
#
# Arguments
#   NONE
#
# Return
#   success - NONE
#   failure - exit 1
#####
RUN_TEST() {
  run_gkpbin $GKP_RUN1_BIN 1
  run_gkpbin $GKP_RUN2_BIN 2

  gen_run1_md5s
  gen_run2_md5s

  if [ ! -f $RUN1_MD5S ]; then
    echo "RUN_TEST: $RUN1_MD5S did not get generated"
    FAILURE
  fi

  if [ ! -f $RUN2_MD5S ]; then
    echo "RUN_TEST: $RUN2_MD5S did not get generated"
    FAILURE
  fi

  ###
  # They should be the same since we are running
  # gatekeeper against the same input used to 
  # generate the baseline MD5s
  ###
  diff -q $RUN1_MD5S $RUN2_MD5S > /dev/null
  ret=$?

  if [ "$ret" -ne "0" ]; then
    ###
    # They were not the same! If the output *should*
    # be different the known good MD5s should be
    # updated!
    ###
    diff $RUN1_MD5S $RUN2_MD5S 
    FAILED
  fi  
}

#####
# MAIN:
#####
if [ "$#" != "6" ]; then
  echo "$usage" 1>&2
  exit 1
fi

###
# Process command line arguments
###
while : ; do
  case "$1" in
    -f )
      shift

      ###
      # If frag file is gbaf.frg then
      # FRAG_PREF is set to gbaf
      ###
      FRAG_PATH=$1
      FRAG_FILE=$(basename $FRAG_PATH)
      ln -s $FRAG_PATH $FRAG_FILE
      FRAG_PREF=${FRAG_FILE%.frg}
      shift     
     ;;
    -a )
      shift
      GKP_RUN1_BIN=$1
      shift
     ;;
    -b )
      shift
      GKP_RUN2_BIN=$1
      shift
     ;;
    * )
      break
     ;;
  esac
done

###
# Check external dependences
###
$(which $XXD_BIN > /dev/null 2>&1) || {
  echo "Err: $XXD_BIN binary not found"
  exit 1
}

###
# INIT
###

###
# BODY
###
RUN_TEST

###
# DONE
###
PASSED

