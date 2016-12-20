gridOptions=-q big.q

genomeSize=4.8m
stopOnReadQuality=false

#-pacbio-raw /data/regression/escherichia_coli_k12.pacbio.p5.fastq.xz

-pacbio-raw /data/reads/escherichia_coli_ne92-p5c3/escherichia_coli_ne92.p5.unknown.m131026_050049_42132_c100598352550000001823109505221424_s1_p0.fasta.xz
-pacbio-raw /data/reads/escherichia_coli_ne92-p5c3/escherichia_coli_ne92.p5.unknown.m131026_082002_42132_c100598352550000001823109505221425_s1_p0.fasta.xz
-pacbio-raw /data/reads/escherichia_coli_ne92-p5c3/escherichia_coli_ne92.p5.unknown.m131026_114301_42132_c100598352550000001823109505221426_s1_p0.fasta.xz
-pacbio-raw /data/reads/escherichia_coli_ne92-p5c3/escherichia_coli_ne92.p5.unknown.m131026_150125_42132_c100598352550000001823109505221427_s1_p0.fasta.xz

onSuccess=/work/canu/src/pipelines/sanity/success.escherichia_coli_ne92.sh

