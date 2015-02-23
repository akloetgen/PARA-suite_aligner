# PARMA
PAr-CLIP Read MApper. Based on the implementation of the Burrows-Wheeler Aligner (BWA; aln/samse/sampe algorithm) Version 0.7.8 and extends it with the use of an empirical error profile.

# Installation
	git clone https://github.com/akloetgen/PARMA.git
	cd PARMA
	make
	
For an easy use of the PARMA algorithm we recommend the additional download of the PARMA toolkit, wrapping the PARMA algorithm into a best-practice pipeline for read mapping (https://github.com/akloetgen/PARMA_tk)

# Usage
	bwa parma [options] <reference_prefix> <input>
	
### Help
	bwa parma
	
### Example on testdata
	bwa index examples/reference_chr1.fa
	bwa parma -p examples/reference_chr1.errorprofile -f examples/testout_simulation_mapped_PARMA.sai examples/reference_chr1.fa examples/testout_simulation.fastq
	bwa samse examples/reference_chr1.fa examples/testout_simulation_mapped_PARMA.sai examples/testout_simulation.fastq > examples/testout_simulation_mapped_PARMA.sam
