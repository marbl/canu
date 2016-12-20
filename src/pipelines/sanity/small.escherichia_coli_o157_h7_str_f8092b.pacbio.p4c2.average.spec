gridOptions=-q big.q

genomeSize=5.4m

stopOnReadQuality=false

#  All files is WAY too much coverage.
#
#  --   Found 1102277 reads.
#  --   Found 4061273543 bases (752.08 times coverage).
#  --
#  --   Read length histogram (one '*' equals 4740.85 reads):
#  --        0    999      0 
#  --     1000   1999 331860 **********************************************************************
#  --     2000   2999 208285 *******************************************
#  --     3000   3999 161792 **********************************
#  --     4000   4999 127542 **************************
#  --     5000   5999  94586 *******************
#  --     6000   6999  66165 *************
#  --     7000   7999  44005 *********
#  --     8000   8999  28529 ******
#  --     9000   9999  17601 ***
#  --    10000  10999  10335 **
#  --    11000  11999   5655 *
#  --    12000  12999   3026 
#  --    13000  13999   1529 
#  --    14000  14999    758 
#  --    15000  15999    370 
#  --    16000  16999    149 
#  --    17000  17999     60 
#  --    18000  18999     21 
#  --    19000  19999      6 
#  --    20000  20999      2 
#  --    21000  21999      1 
#
#  With minReadLength=9999, we get  46.63x coverage, but it fails to generate corrected reads
#
#  With minReadLength=7500, we get 150.40x coverage.
#    Doesn't assemble well, about 30 contigs.  10x in corrected reads.
#
#  With minReadLength=5000, we get 358.04x coverage.
#    Assembles to two contigs, 28x in corrected reads.



#  These are the two largest files, adjacent to each other.  Possibly an unusally good run.
#
#  --   Found 147449 reads.
#  --   Found 613966971 bases (113.69 times coverage).
#  --
#  --   Read length histogram (one '*' equals 511.15 reads):
#  --        0    999      0 
#  --     1000   1999  35781 **********************************************************************
#  --     2000   2999  26363 ***************************************************
#  --     3000   3999  21009 *****************************************
#  --     4000   4999  17417 **********************************
#  --     5000   5999  14247 ***************************
#  --     6000   6999  10635 ********************
#  --     7000   7999   7804 ***************
#  --     8000   8999   5499 **********
#  --     9000   9999   3541 ******
#  --    10000  10999   2336 ****
#  --    11000  11999   1347 **
#  --    12000  12999    747 *
#  --    13000  13999    381 
#  --    14000  14999    184 
#  --    15000  15999    108 
#  --    16000  16999     28 
#  --    17000  17999     16 
#  --    18000  18999      4 
#  --    19000  19999      0 
#  --    20000  20999      2 

#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941219.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941220.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941221.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941222.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941223.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941224.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941225.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941226.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941227.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941228.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941230.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941231.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941233.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941234.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941235.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941242.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941244.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941245.fasta.xz
-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941246.fasta.xz
-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941247.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941248.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941249.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941250.fasta.xz
#-pacbio-raw /data/reads/escherichia_coli_o157:h7_str_f8092b-p4c2/escherichia_coli_o157:h7_str_f8092b.p4c2.nbacc.SRR941251.fasta.xz

onSuccess=/work/canu/src/pipelines/sanity/success.escherichia_coli_o157_h7_str_f8092b.sh
