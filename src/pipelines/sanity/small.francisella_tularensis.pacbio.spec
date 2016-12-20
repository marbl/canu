gridOptions=-q big.q

genomeSize=1.8m

#minReadLength=1000
#stopOnReadQuality=false

#mhapSensitivity=high

#  --   Found 250715 reads.
#  --   Found 741304359 bases (411.83 times coverage).
#  --
#  --   Read length histogram (one '*' equals 1321.45 reads):
#  --        0    999      0 
#  --     1000   1999  92502 **********************************************************************
#  --     2000   2999  63504 ************************************************
#  --     3000   3999  39303 *****************************
#  --     4000   4999  24438 ******************
#  --     5000   5999  14149 **********
#  --     6000   6999   7995 ******
#  --     7000   7999   4387 ***
#  --     8000   8999   2352 *
#  --     9000   9999   1201 
#  --    10000  10999    567 
#  --    11000  11999    197 
#  --    12000  12999     82 
#  --    13000  13999     27 
#  --    14000  14999      8 
#  --    15000  15999      2 
#  --    16000  16999      1 

-pacbio-raw /data/reads/francisella_tularensis/SRR941956_m120810_210002_42144_c100369192550000001523030510101250_s1_p0.fasta.xz
-pacbio-raw /data/reads/francisella_tularensis/SRR941957_m120810_224803_42144_c100369192550000001523030510101251_s1_p0.fasta.xz
-pacbio-raw /data/reads/francisella_tularensis/SRR941958_m120811_003558_42144_c100369192550000001523030510101252_s1_p0.fasta.xz
-pacbio-raw /data/reads/francisella_tularensis/SRR941959_m120811_022401_42144_c100369192550000001523030510101253_s1_p0.fasta.xz
-pacbio-raw /data/reads/francisella_tularensis/SRR941960_m120811_041153_42144_c100369192550000001523030510101254_s1_p0.fasta.xz
-pacbio-raw /data/reads/francisella_tularensis/SRR941961_m120427_215831_42144_c100326282550000001523018509161220_s1_p0.fasta.xz
-pacbio-raw /data/reads/francisella_tularensis/SRR941962_m120428_065000_42144_c100326282550000001523018509161225_s1_p0.fasta.xz
-pacbio-raw /data/reads/francisella_tularensis/SRR941963_m120428_083639_42144_c100326282550000001523018509161226_s1_p0.fasta.xz

onSuccess=/work/canu/src/pipelines/sanity/success.francisella_tularensis.sh
