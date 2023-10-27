#!/bin/bash

#bwa index sacCer3.fa

#for sample in A01_09 A01_11 A01_23 A01_24 A01_27 A01_31 A01_35 A01_39 A01_62 A01_63
#do
#	echo "Aligning sample:" ${sample}
#	bwa mem -t 4 -R "@RG\tID:${sample}\tSM:${sample}" sacCer3.fa ${sample}.fastq > ${sample}.sam

#done


#for sample in A01_09 A01_11 A01_23 A01_24 A01_27 A01_31 A01_35 A01_39 A01_62 A01_63
#do
	#samtools sort -o ${sample}.bam ${sample}.sam
	#samtools index ${sample}.bam
#done

#freebayes -f sacCer3.fa --genotype-qualities -p 1 *.bam > variants.vcf
#vcffilter -f "QUAL > 19" variants.vcf > variants_filtered.vcf

#vcfallelicprimitives -k -g variants_filtered.vcf > variants_decomposed.vcf
#snpEff download R64-1-1.105
#snpEff ann R64-1-1.105 variants_decomposed.vcf > annotated_variants.vcf
#head -n 100 annotated_variants.vcf > shortened_annotation.vcf