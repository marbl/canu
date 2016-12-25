gridOptions=-q big.q

genomeSize=12.1m
stopOnReadQuality=false

-pacbio-raw /data/reads/saccharomyces_cerevisiae_GLBRCY22-3/SRR2989876.fasta.xz
-pacbio-raw /data/reads/saccharomyces_cerevisiae_GLBRCY22-3/SRR2989877.fasta.xz
-pacbio-raw /data/reads/saccharomyces_cerevisiae_GLBRCY22-3/SRR2989878.fasta.xz
-pacbio-raw /data/reads/saccharomyces_cerevisiae_GLBRCY22-3/SRR2989879.fasta.xz
-pacbio-raw /data/reads/saccharomyces_cerevisiae_GLBRCY22-3/SRR2989880.fasta.xz
-pacbio-raw /data/reads/saccharomyces_cerevisiae_GLBRCY22-3/SRR2989881.fasta.xz
-pacbio-raw /data/reads/saccharomyces_cerevisiae_GLBRCY22-3/SRR2989882.fasta.xz
-pacbio-raw /data/reads/saccharomyces_cerevisiae_GLBRCY22-3/SRR2989883.fasta.xz

onSuccess=/work/canu/src/pipelines/sanity/success.saccharomyces_cerevisiae_s288c.sh

