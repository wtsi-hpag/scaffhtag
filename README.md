# scaffhtag v1.0
Pipeline for scaffolding genome assemblies using haplotagging reads.

Pipeline steps:
        
    Scaffolding with scaffhtag:
      1 Barcoded tags are extracted from htag raw sequencing reads and appended 
          to read names for further processing
      2 The reads are mapped to the draft assembly using either BWA or SMALT
      3 Barcodes are sorted together with contigs as well as mapping coordinates
      4 A relation matrix is built to record the shared barcodes among the contigs which may be linked
      5 Order and orientation of linked contigs are determined after nearest neighbours are found. 
      
### Download and Compile:
Requirements for compiling: gcc gcc-4.9.2 or late:

If you see this message,
cc1: error: unrecognised command line option ‘-std=c11’
make: *** [breakhtag.o] Error 1

you need a higher version of gcc
CC= /software/gcc-4.9.2/bin/gcc in the makefile


    $ git clone  https://github.com/wtsi-hpag/scaffhtag.git 
    $ cd htag
    $ bash install.sh
		
If everything compiled successfully you must see the final comment: 
		"Congrats: installation successful!"		


#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) and SMALT (http://www.sanger.ac.uk/science/tools/smalt-0) are downloaded and compiled by htag.

### Run the pipelines

#### Prepare read files with barcode error correction and extraction
           $ /full/path/to/htag/src/scaff_read input.dat htag-reads_RC1.fastq.gz htag-reads_RC2.fastq.gz \
	       input.dat               - input a text file to point the locations of the reads in paired files \
               htag-reads_RC1.fastq.gz - output read file                       \
	       htag-reads_RC1.fastq.gz - output read file                      \

	       input.dat file shoul be like:
		q1=/lustre/scratch116/vr/projects/Tes1_S1_L008_R1_001.fastq.gz \
		q2=/lustre/scratch116/vr/projects/Tes1_S1_L008_R2_001.fastq.gz \
		q1=/lustre/scratch116/vr/projects/Tes1_S2_L008_R1_001.fastq.gz \
		q2=/lustre/scratch116/vr/projects/Tes1_S2_L008_R2_001.fastq.gz \
		q1=/lustre/scratch116/vr/projects/Tes1_S3_L008_R1_001.fastq.gz \
		q2=/lustre/scratch116/vr/projects/Tes1_S3_L008_R2_001.fastq.gz \
		q1=/lustre/scratch116/vr/projects/Tes1_S4_L008_R1_001.fastq.gz \
		q2=/lustre/scratch116/vr/projects/Tes1_S4_L008_R2_001.fastq.gz \
 
#### Run scaffhtag:
           $ /full/path/to/htag/src/scaffhtag -nodes <nodes> -align <aligner> -score <score> \
	   	 -matrix <matrix_size> -read-s1 <min_reads_s1> -read-s2 <min_reads_s2> \
		 -edge <edge_len> -link-s1 <n_links_s1> -link-s2 <n_links_s2> -block <block>  \
		 [ -mkdup Dupmarked.bam ] [ -plot barcode-length.png ] \
		 draft-assembly.fasta htag-reads_RC1.fastq.gz htag-reads_RC2.fastq.gz output_scaffolds.fasta

	       Parameters:
             nodes:        number of CPUs requested  [ default = 30 ]
             score:        averaged mapping score on each barcode fragment [ default = 20 ]
             aligner:      sequence aligner: bwa or smalt [ default = bwa ]
             matrix_size:  relation matrix size [ default = 2000 ]
             min_reads_s1: step 1: minimum number of reads per barcode [ default = 10 ]
             min_reads_s2: step 2: minimum number of reads per barcode [ default = 10 ]
             edge_len:     length of mapped reads to consider for scaffolding [ default = 50000 ]
             n_links_s1:   step 1: minimum number of shared barcodes [ default = 8 ]
             n_links_s2:   step 2: minimum number of shared barcodes [ default = 8 ]
             aggressive:   1 - aggressively mapping filtering on small PacBio/ONT contigs; 
	     		   0 - no aggressive for short read assembly  [ default = 1 ]
             block:        length to determine for nearest neighbours [ default = 50000 ]
             plot:         output image file with barcode length distributions and coverage stats 
	     mkdup:        output bam file with duplicated reads removed \n"); 

	    


