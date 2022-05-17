# scaffhtag v2.0
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
           $ /full/path/to/htag/src/scaff_read input.dat htag-reads_BC1.fastq.gz htag-reads_BC2.fastq.gz \
	       input.dat               - input a text file to point the locations of the reads in cram files \
               htag-reads_BC1.fastq.gz - output read file                       \
	       htag-reads_BC1.fastq.gz - output read file                      \

	       input.dat file shoul be like with full path:
		/lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/43969#17.cram \
		/lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/43969#18.cram \
		/lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/43969#19.cram \
		/lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/43969#20.cram \
		/lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/43969#21.cram \
		/lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/43969#22.cram \
		/lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/43969#23.cram \
		/lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/43969#24.cram \ 

#### Run scaffhtag with paired reads htag-reads_BC1.fastq.gz htag-reads_BC2.fastq.gz:
           $ /full/path/to/htag/src/scaffhtag -nodes <nodes> -align <aligner> -score <score> \
	   	 -matrix <matrix_size> -read-s1 <min_reads_s1> -read-s2 <min_reads_s2> \
		 -edge <edge_len> -link-s1 <n_links_s1> -link-s2 <n_links_s2> -block <block>  \
		 [ -mkdup Dupmarked.bam ] [ -plot barcode-length.png ] \
		 draft-assembly.fasta htag-reads_BC1.fastq.gz htag-reads_BC2.fastq.gz output_scaffolds.fasta

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


#### Run scaffhtag with aligned and sorted bam file: aligned.bam  
           $ /full/path/to/htag/src/scaffhtag -nodes <nodes> -plot barcode-length.png -bam /lustre/scratch117/sciops/team117/hpag/zn1/aligned.bam draft-assembly.fasta output_scaffolds.fasta \


#### Run alignment with ema:
           $ /full/path/to/htag/src/scaff-bin/ema-align.csh <input_cram_file> <Output_workdirectory> <bwa_index> <output_bam_file> \
	   
#####	 Instructions for Installation
 	   https://github.com/wtsi-hpag/scaffhtag
   	 chmod 777 ema-align.csh 
	   Tools needed 
	     1. samtools version 1.15 or later 
	     2. SamHaplotag 
 	        https://github.com/wtsi-hpag/SamHaplotag  
	     3. bwa version 0.7.12-r1044 or later
	     4. ema 
	        https://github.com/arshajii/ema 
	 Use of bioconda for installation 
	 Reference index 
	 -- Say you have a reference genome assembly  
	   /lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/bindex/Oak-chr.fasta
	   cd /lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/bindex/
	   bwa index Oak-chr.fasta
	   samtools faidx Oak-chr.fasta
	 Alignment 
	 -- Say you have cram file 43969#17.cram and Oak-chr.fasta index 
	   /nfs/users/nfs_z/zn1/bin/ema-align.csh 43969#17.cram readsplit-17 /lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/bindex/Oak-chr.fasta ema_final-17.bam
	 Your output file ema_final-17.bam will be in readsplit-17.  
  
 


