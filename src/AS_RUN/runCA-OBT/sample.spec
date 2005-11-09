
############################################################
#
#  Job Parameters
#
useGrid = 0


#  The location of your wgs-assembler directory, e.g.,
#  $binRoot/Linux/bin/consensus
#
#  If not set, it defaults to the location where the runCA-OBT.pl
#  script is.
#
#binRoot = /bioinfo/work/projects/macaque-v2/wgs


#  Host specification.  The local host is used for large memory and
#  sequential tasks, the grid host is used for parallel steps.
#
#  Valid options here are
#
#    Machine   Host
#    i686      Linux       #  32 bit linux
#    x86_64    Linux       #  64 bit linux on opteron
#    i386      FreeBSD     #  32 bit FreeBSD
#    --        OSF1        #  64 bit compaq OSF1 on alpha
#    alpha     --
#    --        Darwin      #  32/64 bit OS-X (on PPC)
#
localMachine = x86_64
localHost    = Linux
gridMachine  = i686
gridHost     = Linux


############################################################
#
#  Overlap based trim parameters
#
doOverlapTrimming = 1

#  OBT works best if it knows where the vector is.
#
#vectorIntersect = /path/to/file

#  OBT needs to sort all the overlaps, more memory here is better.
#  Probably should be the same as ovlStoreMemory.
#
ovlSortMemory = 1000


############################################################
#
#  Overlapper parameters
#
#  ovlThreads -- number of threads per overlapper job.  On the VI
#  grid, hosts are hyper-threaded dual-Xeons.  Previous trivial tests
#  indicated that using two threads per CPU gave about 20% better
#  performance than one thread per CPU.
#
#  ovlHashBlockSize -- number of fragments to build the in-core hash
#  table with.  200,000 fragments of mean size ~550bp ran in 4GB at
#  Celera.  VI fragments are ~800bp.  A host at VI has 2GB, and can
#  run two processes, so figure on ~600MB per job.  The script for
#  running dog used 40000.  30000 was used here to prevent any chance
#  of paging.
#
#  ovlRefBlockSize -- to better utilize CPU and to make jobs shorter, we
#  can segment the number of fragments we run by the overlapper in one
#  run.  1,000,000 (with 40,000 ovlHaskBlockSize) is reported to give
#  about one hour of run time.
#
ovlThreads       = 2
ovlHashBlockSize = 40000
ovlRefBlockSize  = 2000000
ovlMemory        = 1GB

ovlThreads       = 2
ovlHashBlockSize = 150000
ovlRefBlockSize  = 5000000
ovlMemory        = 2GB

#  Amount of memory to use when building the overlap store.  More is
#  always better.  Probably should be the same as ovlStoreMemory.
#
ovlStoreMemory   = 1000

############################################################
#
#  Fragment correction parameters
#
#  Using the store directly (recent commit to the tree, compile time
#  option), testing shows that, for a human assembly using 22-mer
#  and standard overlaps:
#
#      10,000 frags per batch needs   132 MB
#      50,000 frags per batch needs   650 MB
#     100,000 frags per batch needs  1300 MB
#     200,000 frags per batch needs  2500 MB
#     500,000 frags per batch needs  6300 MB
#   1,000,000 frags per batch needs 13000 MB
#   2,000,000 frags per batch needs 23000 MB
#   2,500,000 frags per batch needs 30000 MB (died)
#
#  3 million fragments work on a 32GB box (usually), assuming
#  correct-frags doesn't use a temporary internal store.
#
frgCorrBatchSize = 175000
frgCorrThreads   = 2
frgCorrOnGrid    = 0


############################################################
#
#  Overlap correction parameters
#
#  Don't know anything about sizes here....
#
ovlCorrBatchSize = 175000
ovlCorrOnGrid    = 0


############################################################
#
#  Unitigger parameters
#
#  The default values are usually acceptable.
#
#genomeSize         = -l 185000000
#unitiggerEdges     = -m 95000000
#unitiggerFragments = -n 30000000


############################################################
#
#  Scaffolding parameters
#

#  The stone level of the final cgw invocation.
#
stoneLevel       = 2

#  If enabled, it will do $eCRRounds iterations of extendClearRanges.
#  eCRRounds should be at least one.
#
doExtendClearRanges = 0
eCRRounds = 1
