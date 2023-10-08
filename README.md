# master_thesis_2023
Load the packages:
- conda activate snakemake
- module load FASTQC
- module load SAMtools
1. Run Snakefile:
- snakemake --cores 24 --j 24
- Note: fastqc and multiqc - 1 thread per 24 jobs; bwa-mem2 - 24 threads per 1 job
