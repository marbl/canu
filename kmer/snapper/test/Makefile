


all:
	@../../fasta/leaff -G 1 2000 2000 > 1.f
	@../../fasta/leaff -G 1  700  700 > 2.f
	@../../fasta/leaff -G 1 3000 3000 > 3.f
	@echo ">g1" > g.f
	@cat 1.f 2.f 3.f | grep -v '>' >> g.f
	@echo ">g2" >> g.f
	@cat 3.f 2.f 1.f | grep -v '>' >> g.f
	@echo ">g3" >> g.f
	@cat 1.f 2.f 3.f | grep -v '>' >> g.f
	@rm -f 1.f 3.f
	../snapper2 -mersize 23 -queries 2.f -genomic g.f -tmp x.tmp -output x.out
	@cat x.out
	@echo ""
	@echo "Check that the ranges are 1900-2800, 2900-3800 and 1900-2800"
	@echo "Check that ddefs change"
	@echo "Check that exons are the same"
	@rm -f 1.f* 2.f* 3.f* g.f* x.tmp* x.out*
