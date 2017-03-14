gridOptions=-q big.q

genomeSize=100m

minReadLength=12500        #  Gives 55.88x
stopOnReadQuality=false

#
#  -- In gatekeeper store 'correction/test.gkpStore':
#  --   Found 700739 reads.
#  --   Found 8095619156 bases (80.96 times coverage).
#  --
#  --   Read length histogram (one '*' equals 667.82 reads):
#  --        0    999      0 
#  --     1000   1999  46748 **********************************************************************
#  --     2000   2999  42946 ****************************************************************
#  --     3000   3999  39761 ***********************************************************
#  --     4000   4999  37550 ********************************************************
#  --     5000   5999  35871 *****************************************************
#  --     6000   6999  34847 ****************************************************
#  --     7000   7999  34368 ***************************************************
#  --     8000   8999  32753 *************************************************
#  --     9000   9999  30850 **********************************************
#  --    10000  10999  28637 ******************************************
#  --    11000  11999  27052 ****************************************
#  --    12000  12999  25478 **************************************
#  --    13000  13999  24504 ************************************
#  --    14000  14999  24549 ************************************
#  --    15000  15999  26035 **************************************
#  --    16000  16999  31150 **********************************************
#  --    17000  17999  33713 **************************************************
#  --    18000  18999  29128 *******************************************
#  --    19000  19999  23535 ***********************************
#  --    20000  20999  18692 ***************************
#  --    21000  21999  14486 *********************
#  --    22000  22999  11399 *****************
#  --    23000  23999   9051 *************
#  --    24000  24999   7273 **********
#  --    25000  25999   5743 ********
#  --    26000  26999   4722 *******
#  --    27000  27999   3883 *****
#  --    28000  28999   2963 ****
#  --    29000  29999   2455 ***
#  --    30000  30999   2087 ***
#  --    31000  31999   1604 **
#  --    32000  32999   1384 **
#  --    33000  33999   1207 *
#  --    34000  34999    976 *
#  --    35000  35999    766 *
#  --    36000  36999    612 
#  --    37000  37999    465 
#  --    38000  38999    401 
#  --    39000  39999    306 
#  --    40000  40999    228 
#  --    41000  41999    143 
#  --    42000  42999    119 
#  --    43000  43999     74 
#  --    44000  44999     69 
#  --    45000  45999     53 
#  --    46000  46999     32 
#  --    47000  47999     20 
#  --    48000  48999     16 
#  --    49000  49999      4 
#  --    50000  50999      9 
#  --    51000  51999     12 
#  --    52000  52999      5 
#  --    53000  53999      2 
#  --    54000  54999      2 
#  --    55000  55999      1 
#

-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_062819_ethan_c100699582550000001823139903261540_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_062819_ethan_c100699582550000001823139903261540_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_062819_ethan_c100699582550000001823139903261540_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_104939_ethan_c100699582550000001823139903261541_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_104939_ethan_c100699582550000001823139903261541_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_104939_ethan_c100699582550000001823139903261541_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_151111_ethan_c100699582550000001823139903261542_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_151111_ethan_c100699582550000001823139903261542_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_151111_ethan_c100699582550000001823139903261542_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_184123_42139_c100719602550000001823155305141590_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_184123_42139_c100719602550000001823155305141590_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_184123_42139_c100719602550000001823155305141590_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_191128_sidney_c100699772550000001823139903261590_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_191128_sidney_c100699772550000001823139903261590_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_191128_sidney_c100699772550000001823139903261590_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_192713_ethan_c100699582550000001823139903261543_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_192713_ethan_c100699582550000001823139903261543_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_192713_ethan_c100699582550000001823139903261543_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_230547_42139_c100719602550000001823155305141591_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_230547_42139_c100719602550000001823155305141591_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_230547_42139_c100719602550000001823155305141591_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_233028_sidney_c100699772550000001823139903261591_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_233028_sidney_c100699772550000001823139903261591_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_233028_sidney_c100699772550000001823139903261591_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_234420_ethan_c100699582550000001823139903261544_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_234420_ethan_c100699582550000001823139903261544_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140928_234420_ethan_c100699582550000001823139903261544_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_033247_42139_c100719602550000001823155305141592_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_033247_42139_c100719602550000001823155305141592_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_033247_42139_c100719602550000001823155305141592_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_034941_sidney_c100699772550000001823139903261592_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_034941_sidney_c100699772550000001823139903261592_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_034941_sidney_c100699772550000001823139903261592_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_040333_ethan_c100699582550000001823139903261545_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_040333_ethan_c100699582550000001823139903261545_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_040333_ethan_c100699582550000001823139903261545_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_075857_42139_c100719602550000001823155305141593_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_075857_42139_c100719602550000001823155305141593_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_075857_42139_c100719602550000001823155305141593_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_080908_sidney_c100699772550000001823139903261593_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_080908_sidney_c100699772550000001823139903261593_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_080908_sidney_c100699772550000001823139903261593_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_082928_ethan_c100699582550000001823139903261546_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_082928_ethan_c100699582550000001823139903261546_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_082928_ethan_c100699582550000001823139903261546_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_122654_42139_c100719602550000001823155305141594_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_122654_42139_c100719602550000001823155305141594_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_122654_42139_c100719602550000001823155305141594_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_122826_sidney_c100699772550000001823139903261594_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_122826_sidney_c100699772550000001823139903261594_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_122826_sidney_c100699772550000001823139903261594_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_124838_ethan_c100699582550000001823139903261547_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_124838_ethan_c100699582550000001823139903261547_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_124838_ethan_c100699582550000001823139903261547_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_164720_sidney_c100699772550000001823139903261595_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_164720_sidney_c100699772550000001823139903261595_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_164720_sidney_c100699772550000001823139903261595_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_211052_sidney_c100699772550000001823139903261596_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_211052_sidney_c100699772550000001823139903261596_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140929_211052_sidney_c100699772550000001823139903261596_s1_p0.3.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140930_013011_sidney_c100699772550000001823139903261597_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140930_013011_sidney_c100699772550000001823139903261597_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/regression/reads/caenorhabditis_elegans-p6c4/m140930_013011_sidney_c100699772550000001823139903261597_s1_p0.3.subreads.fasta.xz

onSuccess=/work/canu/src/pipelines/sanity/success.caenorhabditis_elegans.sh
