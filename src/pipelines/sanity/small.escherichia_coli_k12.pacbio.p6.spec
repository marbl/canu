gridOptions=-q big.q

genomeSize=4.8m
stopOnReadQuality=false

#-pacbio-raw /data/regression/escherichia_coli_k12.pacbio.p6.fastq.xz

-pacbio-raw /data/reads/escherichia_coli_k12/escherichia_coli_k12_p6c4.m140930_121059_sherri_c100688052550000001823139503241542_s1_p0.1.subreads.fasta.xz
-pacbio-raw /data/reads/escherichia_coli_k12/escherichia_coli_k12_p6c4.m140930_121059_sherri_c100688052550000001823139503241542_s1_p0.2.subreads.fasta.xz
-pacbio-raw /data/reads/escherichia_coli_k12/escherichia_coli_k12_p6c4.m140930_121059_sherri_c100688052550000001823139503241542_s1_p0.3.subreads.fasta.xz

onSuccess=/work/canu/src/pipelines/sanity/success.escherichia_coli_k12.sh

