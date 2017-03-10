gridOptions=-q big.q

genomeSize=4.8m
stopOnReadQuality=false

#
#  --   Found 7937 reads.
#  --   Found 1142222431 bases (237.96 times coverage).
#  --
#  --   Read length histogram (one '*' equals 2.75 reads):
#
#  Min read length = 100000
#  Max read length = 790000
#

-nanopore-raw /data/regression/reads/escherichia_coli_k12.nanopore.r9.4.superlong.fasta.xz

onSuccess=/work/canu/src/pipelines/sanity/success.escherichia_coli_k12.sh

