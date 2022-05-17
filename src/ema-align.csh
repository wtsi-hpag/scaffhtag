#!/bin/tcsh
# by zn1@sanger.ac.uk
# 05/05/2022


# Instructions for Installation
# https://github.com/wtsi-hpag/scaffhtag
# chmod 777 ema-align.csh 
# Tools needed 
#   1. samtools version 1.15 or later 
#   2. SamHaplotag 
#      https://github.com/wtsi-hpag/SamHaplotag  
#   3. bwa version 0.7.12-r1044 or later
#   4. ema 
#      https://github.com/arshajii/ema 
# Use of bioconda for installation 
# Reference index 
# -- Say you have a reference genome assembly  
#   /lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/bindex/Oak-chr.fasta
#   cd /lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/bindex/
#   bwa index Oak-chr.fasta
#   samtools faidx Oak-chr.fasta
# Alignment 
# -- Say you have cram file 43969#17.cram and Oak-chr.fasta index 
#   /nfs/users/nfs_z/zn1/bin/ema-align.csh 43969#17.cram newsplit-17 /lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/bindex/Oak-chr.fasta ema_final-17.bam
# Your output file ema_final-17.bam will be in newsplit-17.  
#  


if($#argv < 2) then
        echo "Usage:    $0 [options] cram_file work_directory refer_index output_file"
        exit
endif

set cram_file=$1
set work_dirt=$2
set refe_dirt=$3
set outp_file=$4

samtools view -h@ 16 -F 0xF00 $cram_file | SamHaplotag | samtools view -@ 16 -o tagged_reads.cram
samtools fastq -@ 16 -nT BX -0 /dev/null -s /dev/null tagged_reads.cram | 10xSpoof SamHaplotag_Clear_BC | bgzip -@ 16 >10x_spoofed_reads.fq.gz

samtools fastq -@ 16 -nT BX -0 /dev/null -s /dev/null tagged_reads.cram | 16BaseBCGen | bgzip -@ 16 >16BaseBC_reads.fq.gz
cut -f 2 HaploTag_to_16BaseBCs | tail -n +2 >16BaseBCs


zcat 16BaseBC_reads.fq.gz | ema count -w 16BaseBCs -o new-readcount
zcat 16BaseBC_reads.fq.gz | ema preproc -w 16BaseBCs -n 1 -o $work_dirt new-readcount.ema-ncnt


cd $work_dirt 

#/lustre/scratch117/sciops/team117/hpag/zn1/project/HiC/QC/run-43969/oak1/bindex/Oak-chr.fasta
bwa mem -p -t 50 -M -R "@RG\tID:rg1\tSM:sample1" $refe_dirt ema-nobc > ema-nobc.sam

samtools view -@ 50 -bS ema-nobc.sam > ema-nobc.bam 
samtools sort -@ 50 -o ema-nobc-sort.bam -O BAM ema-nobc.bam > try.out

ema align -t 50 -d -r $refe_dirt -o align.sam -s ema-bin-000 > try.out

samtools view -bS align.sam > align.bam
samtools sort -@ 50 -o align-sorted.bam -O BAM align.bam > try.out
rm -rf align.sam align.bam ema-nobc.sam ema-nobc.bam 

samtools merge -@ 50 -o $outp_file align-sorted.bam ema-nobc-sort.bam 


