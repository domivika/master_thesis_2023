# master_thesis_2023
Load the packages:
- conda activate thesis #snakemake; sambamba; multiqc; fastqc; picard; gatk; samtools
- module load picard
- 
1. Run Snakefile for single-cell Illumina:
- snakemake --cores 24 -j 24
- Note: fastqc and multiqc - 1 thread per 24 jobs; bwa-mem2 - 24 threads per 1 job

2. Run Snakefile for bulk Illumina:
- snakemake --cores 24 -j 1
- Note: fastqc, multiqc and bwa-mem2 - 24 threads per 1 job

#TODOs:
- Fix FastQC rule to  run on all samples in parallel (currently it's running one sample at a time)
- Define the output of GenomicsDBImport
- Run HaplotypeCaller in parallel on multiple chromosomes
- Add log files for each rule
- Check MergeSortedBAM if can be writen in a prettier way
- Check all the rules if they can be more parallized  
