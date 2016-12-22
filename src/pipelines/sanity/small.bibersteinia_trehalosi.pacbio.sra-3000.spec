gridOptions=-q big.q

genomeSize=2.5m

minReadLength=3000
stopOnReadQuality=false

mhapSensitivity=high

#  --   Found 269835 reads.
#  --   Found 1310159772 bases (524.06 times coverage).
#  --
#  --   Read length histogram (one '*' equals 1453.42 reads):
#  --        0    999      0 
#  --     1000   1999      0 
#  --     2000   2999      0 
#  --     3000   3999 101740 **********************************************************************
#  --     4000   4999  70730 ************************************************
#  --     5000   5999  43872 ******************************
#  --     6000   6999  25286 *****************
#  --     7000   7999  13648 *********
#  --     8000   8999   7356 *****
#  --     9000   9999   3671 **
#  --    10000  10999   1820 *
#  --    11000  11999    843 
#  --    12000  12999    468 
#  --    13000  13999    238 
#  --    14000  14999    102 
#  --    15000  15999     45 
#  --    16000  16999     12 
#  --    17000  17999      3 
#  --    18000  18999      1 

-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849064.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849065.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849066.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849067.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849068.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849069.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849073.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849074.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849075.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849076.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849077.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849078.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849079.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849080.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849081.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849082.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849083.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849084.fastq.xz
-pacbio-raw /data/reads/bibersteinia_trehalosi_USDA-ARS-USMARC-192-p4c2/SRR849085.fastq.xz

onSuccess=/work/canu/src/pipelines/sanity/success.bibersteinia_trehalosi.sh

