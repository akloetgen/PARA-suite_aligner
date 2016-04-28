# PARA-suite aligner
PAR-CLIP Analyzing suite. Based on the implementation of the Burrows-Wheeler Aligner (BWA; aln/samse/sampe algorithm) Version 0.7.8 and extends it with the use of an empirical error profile.

# Installation
	git clone https://github.com/akloetgen/PARA-suite_aligner.git
	cd PARA-suite_aligner
	make
	
For an easy use of the PARA-suite ailgner we recommend the additional download of the PARA-suite, wrapping the PARA-suite aligner into a best-practice pipeline for read mapping (https://github.com/akloetgen/PARA-suite)

# Usage
	./bwa parasuite [options] <reference_fasta> <input>
	
### Help
	./bwa parasuite
	
### Example on testdata
	./bwa index examples/reference_chr1.fa
	./bwa parasuite -p examples/reference_chr1.errorprofile -f examples/testout_simulation_mapped_PARMA.sai examples/reference_chr1.fa examples/testout_simulation.fastq
	./bwa samse examples/reference_chr1.fa examples/testout_simulation_mapped_PARMA.sai examples/testout_simulation.fastq > examples/testout_simulation_mapped_PARMA.sam
