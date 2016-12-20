gridOptions=-q big.q

genomeSize=4.8m
stopOnReadQuality=false

#-pacbio-raw /data/regression/escherichia_coli_k12.pacbio.p4.fastq.xz

-pacbio-raw /data/reads/escherichia_coli_ne92-p4c2/escherichia_coli_ne92.p4.unknown.m131113_220250_42132_c100587992550000001823102004281411_s1_p0.fasta.xz
-pacbio-raw /data/reads/escherichia_coli_ne92-p4c2/escherichia_coli_ne92.p4.unknown.m131114_002159_42132_c100587992550000001823102004281412_s1_p0.fasta.xz

onSuccess=/work/canu/src/pipelines/sanity/success.escherichia_coli_ne92.sh

