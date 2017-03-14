gridOptions=-q big.q

genomeSize=4.8m
stopOnReadQuality=false

#
#  -- In gatekeeper store 'correction/small.escherichia_coli_ne92.pacbio.p5.gkpStore':
#  --   Found 401653 reads.
#  --   Found 2279072149 bases (474.8 times coverage).
#  --
#  --   Read length histogram (one '*' equals 996.88 reads):
#  --        0    999      0 
#  --     1000   1999  69782 **********************************************************************
#  --     2000   2999  64880 *****************************************************************
#  --     3000   3999  52563 ****************************************************
#  --     4000   4999  38134 **************************************
#  --     5000   5999  30848 ******************************
#  --     6000   6999  25500 *************************
#  --     7000   7999  21151 *********************
#  --     8000   8999  18286 ******************
#  --     9000   9999  15887 ***************
#  --    10000  10999  13965 **************
#  --    11000  11999  11796 ***********
#  --    12000  12999   9706 *********
#  --    13000  13999   7986 ********
#  --    14000  14999   6084 ******
#  --    15000  15999   4545 ****
#  --    16000  16999   3432 ***
#  --    17000  17999   2441 **
#  --    18000  18999   1770 *
#  --    19000  19999   1137 *
#  --    20000  20999    747 
#  --    21000  21999    444 
#  --    22000  22999    263 
#  --    23000  23999    149 
#  --    24000  24999     86 
#  --    25000  25999     37 
#  --    26000  26999     18 
#  --    27000  27999     10 
#  --    28000  28999      3 
#  --    29000  29999      1 
#  --    30000  30999      2 
#

-pacbio-raw /data/regression/reads/escherichia_coli_ne92-p5c3/escherichia_coli_ne92.p5.unknown.m131026_050049_42132_c100598352550000001823109505221424_s1_p0.fasta.xz
-pacbio-raw /data/regression/reads/escherichia_coli_ne92-p5c3/escherichia_coli_ne92.p5.unknown.m131026_082002_42132_c100598352550000001823109505221425_s1_p0.fasta.xz
-pacbio-raw /data/regression/reads/escherichia_coli_ne92-p5c3/escherichia_coli_ne92.p5.unknown.m131026_114301_42132_c100598352550000001823109505221426_s1_p0.fasta.xz
-pacbio-raw /data/regression/reads/escherichia_coli_ne92-p5c3/escherichia_coli_ne92.p5.unknown.m131026_150125_42132_c100598352550000001823109505221427_s1_p0.fasta.xz

onSuccess=/work/canu/src/pipelines/sanity/success.escherichia_coli_ne92.sh

