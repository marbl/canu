
passed: ../positionDB ../../fasta/leaff
	../../fasta/leaff -G 1000 5000 10000 > x.fasta
	../positionDB -test1 x.fasta
	../../fasta/leaff -G 100 500 1000 > y.fasta
	../positionDB -test2 x.fasta y.fasta
	rm -f x.fasta y.fasta
	touch passed

clean:
	rm -f passed x.fasta y.fasta
