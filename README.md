# PARMA
PAr-CLIP Read MApper. Based on the implementation of the Burrows-Wheeler Aligner (BWA; aln/samse/sampe algorithm) Version 0.7.8 and extends it with the use of an empirical error profile.

# Installation
	git clone https://github.com/akloetgen/PARMA.git
	cd PARMA
	make
	
For an easy use of the PARMA algorithm we recommend the additional download of the PARMA toolkit, wrapping the PARMA algorithm into a best-practice pipeline for read mapping (https://github.com/akloetgen/PARMA_tk)

# Usage
	bwa parma [options] <reference_prefix> <input>
	
# Example on testdata
	bwa parma -p examples/references/reference_chr1.errorprofile -f examples/mapping/testout_simulation_mapped_PARMA.sai examples/references/reference_chr1.fa examples/simulation/testout_simulation.fastq
	bwa samse examples/references/reference_chr1.fa examples/mapping/testout_simulation_mapped_PARMA.sai examples/simulation/testout_simulation.fastq > examples/mapping/testout_simulation_mapped_PARMA.sam
