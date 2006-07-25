########################################
#
#  See doc.tex for descriptions.
#


############################################################
#
#  Job Parameters
#
useGrid   = 0
binRoot   = undef
scratch   = /scratch

#localMachine = x86_64
#localHost    = Linux
#gridMachine  = i686
#gridHost     = Linux

fakeUIDs      = 0
uidServer     = undef


############################################################
#
#  Overlap based trim parameters
#
doOverlapTrimming = 1
doBackupFragStore = 1
vectorIntersect   = undef
immutableFrags    = undef
ovlSortMemory     = 1024


############################################################
#
#  Overlapper parameters
#
ovlOnGrid        = 1
ovlThreads       = 2
ovlHashBlockSize = 40000
ovlRefBlockSize  = 2000000
ovlMemory        = 1GB
ovlStoreMemory   = 1024


############################################################
#
#  Fragment correction parameters
#
doFragmentCorrection = 1
frgCorrOnGrid        = 0
frgCorrBatchSize     = 175000
frgCorrThreads       = 2
frgCorrConcurrency   = 1
ovlCorrOnGrid        = 0
ovlCorrBatchSize     = 175000
ovlCorrConcurrency   = 4


############################################################
#
#  Unitigger parameters
#
#  The default values are usually acceptable.
#
utgEdges          = undef
utgErrorRate      = 15
utgFragments      = undef
utgBubblePopping  = 1
utgGenomeSize     = undef


############################################################
#
#  Consensus parameters
#
cnsPartitions    = 128
cnsMinFrags      = 75000
cnsConcurrency   = 2
cnsOnGrid        = 1


############################################################
#
#  Scaffolding parameters
#
delayInterleavedMerging = 0
stoneLevel              = 2
doExtendClearRanges     = 2
updateDistanceType      = "pre"
doUpdateDistanceRecords = 1
doResolveSurrogates     = 1
