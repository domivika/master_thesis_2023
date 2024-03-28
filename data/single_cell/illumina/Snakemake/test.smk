# Define paths to references and directories
reference = "/lustre1/project/stg_00096/references/GRCh38.alt-masked-V2/index/bwa-mem2/Homo_sapiens_assembly38_masked.fasta"
known_sites = "/staging/leuven/stg_00096/references/GRCh38.alt-masked-V2/annotation/snv_indel/gnomAD_v3.1.2/small-gnomad-common-GRCexcl.genomes.v3.1.2.sites.all.vcf.bgz"
GATK_path = "/staging/leuven/stg_00096/software/singularity/gatk_4.2.6.1.sif"
work_dir = "/staging/leuven/stg_00096/home/dkresa/ResolveOME/single_cell/230623.NovaSeq2.FCA/1277"
stage_dir = "/staging/leuven/stg_00096"
lustre_dir = "/lustre1/project/stg_00096"
scratch_dir = "/scratch/leuven/343/vsc34319"

# Define global wildcards
files = glob_wildcards("/staging/leuven/stg_00096/home/dkresa/ResolveOME/single_cell/230623.NovaSeq2.FCA/1277/{sample}_L00{l}_R{r}_001.fastq.gz")
fq1 = glob_wildcards("/staging/leuven/stg_00096/home/dkresa/ResolveOME/single_cell/230623.NovaSeq2.FCA/1277/{sample}_L00{l}_R1_001.fastq.gz")
fq2 = glob_wildcards("/staging/leuven/stg_00096/home/dkresa/ResolveOME/single_cell/230623.NovaSeq2.FCA/1277/{sample}_L00{l}_R2_001.fastq.gz")
chrom = glob_wildcards("/staging/leuven/stg_00096/home/dkresa/data/single_cell/230623.NovaSeq2.FCA/Snakemake/contigs/chr{n}")


# Define output files
fastqc_output = expand("QC/1_premap/fastqc/{sample}_L00{l}_R{r}_001_fastqc.html", zip, sample=files.sample, l=files.l, r=files.r)
multiqc_output = "QC/1_premap/multiqc/multiqc_report.html"
bwa_mem2_output = expand("raw_bams/{sample}_L00{l}.bam", zip, sample=files.sample, l=files.l)
sort_bam_files_output = expand("raw_bams/{sample}_L00{l}.sorted.bam", zip, sample=files.sample, l=files.l)
mark_duplicates_output = expand("raw_bams/{sample}_L00{l}.sorted.marked.bam", zip, sample=files.sample, l=files.l)
base_recalibrator_output = expand("recalibrated/{sample}_L00{l}.recal_data.table", zip, sample=files.sample, l=files.l)
apply_bqsr_output = expand("recalibrated/{sample}_L00{l}.sorted.marked.recal.bam", zip, sample=files.sample, l=files.l)
merge_sorted_bam_output = expand("bam/{sample}.bam", sample=files.sample)
haplotype_calling_output = expand("gvcf/{sample}/{sample}_chr{n}.g.vcf.gz", sample=files.sample, n=chrom.n)
GatherVcfs_output = expand("GatherVcfs/{sample}.g.vcf.gz", sample=files.sample)
GenomicsDBImport_output = expand("database/chr{n}_gdb", n=chrom.n)
GenotypeGVCFs_output = expand("vcf/chr{n}.vcf", n=chrom.n)
GatherVcfs2_output = "vcf/merged.vcf"
VariantFiltration_output = "vcf/merged_filtered.vcf"
CollectWgsMetrics_output = expand("QC/2_postmap/{sample}_wgs_metrics.txt", sample=files.sample)
CollectAlignmentSummaryMetrics_output = expand("QC/2_postmap/{sample}_alignment_summary_metrics.txt", sample=files.sample)

# Define output reports
bwa_mem2_report = expand("raw_bams/reports/bwa_mem2/{sample}_L00{l}.txt", zip, sample=files.sample, l=files.l)
sort_bam_files_report = expand("raw_bams/reports/sort_bam_files/{sample}_L00{l}.txt", zip, sample=files.sample, l=files.l)
mark_duplicates_report = expand("raw_bams/reports/mark_duplicates/{sample}_L00{l}.txt", zip, sample=files.sample, l=files.l)
base_recalibrator_report = expand("recalibrated/reports/base_recalibrator/{sample}_L00{l}.txt", zip, sample=files.sample, l=files.l)
apply_bqsr_report = expand("recalibrated/reports/apply_bqsr/{sample}_L00{l}.txt", zip, sample=files.sample, l=files.l)
merge_sorted_bam_report = expand("bam/reports/merge_sorted_bam/{sample}.txt", sample=files.sample)
haplotype_calling_report = expand("reports/haplotype_calling/{sample}/{sample}_chr{n}.txt", sample=files.sample, n=chrom.n)
GatherVcfs_report = expand("reports/GatherVcfs/{sample}.txt", sample=files.sample)
GenomicsDBImport_report = expand("reports/GenomicsDBImport/chr{n}.txt", n=chrom.n)
GenotypeGVCFs_report = expand("reports/GenotypeGVCFs/chr{n}.txt", n=chrom.n)
GatherVcfs2_report = "reports/GatherVcfs2/merged.txt"
VariantFiltration_report = "reports/VariantFiltration/merged_filtered.txt"
CollectWgsMetrics_report = expand("reports/CollectWgsMetrics/{sample}_wgs_metrics.txt", sample=files.sample)
CollectAlignmentSummaryMetrics_report = expand("reports/CollectAlignmentSummaryMetrics/{sample}_alignment_summary_metrics.txt", sample=files.sample)


rule all:
    input:
        fastqc_output,
        multiqc_output,
        bwa_mem2_output,
        bwa_mem2_report,
        sort_bam_files_output,
        sort_bam_files_report,
        mark_duplicates_output,
        mark_duplicates_report,
        base_recalibrator_output,
        base_recalibrator_report,
        apply_bqsr_output,
        apply_bqsr_report,
        merge_sorted_bam_output,
        merge_sorted_bam_report,
        haplotype_calling_output,
        haplotype_calling_report,
        GatherVcfs_output,
        GatherVcfs_report,
        GenomicsDBImport_output,
        GenomicsDBImport_report,
        GenotypeGVCFs_output,
        GenotypeGVCFs_report,
        GatherVcfs2_output,
        GatherVcfs2_report,
        VariantFiltration_output,
        VariantFiltration_report,
        CollectWgsMetrics_output,
        CollectWgsMetrics_report,
        CollectAlignmentSummaryMetrics_output,
        CollectAlignmentSummaryMetrics_report

# Rule for quality control (FastQC):
rule fastqc:
    input:
        fq = "/staging/leuven/stg_00096/home/dkresa/ResolveOME/single_cell/230623.NovaSeq2.FCA/1277/{sample}_L00{l}_R{r}_001.fastq.gz"
    output:
        zip = "QC/1_premap/fastqc/{sample}_L00{l}_R{r}_001_fastqc.zip",
        html = "QC/1_premap/fastqc/{sample}_L00{l}_R{r}_001_fastqc.html"
    threads: 1
    shell:
        """
        echo "Input Fastq: {input.fq}"
        fastqc -o QC/1_premap/fastqc {input.fq}
        """

# Rule for generating MultiQC report
rule multiqc:
    input:
        report_files = fastqc_output
    output:
        html_report = "QC/1_premap/multiqc/multiqc_report.html"
    threads: 24
    shell:
        """
        multiqc QC/1_premap/fastqc -o QC/1_premap/multiqc
        """

# Rule for BWA MEM2 alignment
rule bwa_mem2:
    input:
        fq1 = "/staging/leuven/stg_00096/home/dkresa/ResolveOME/single_cell/230623.NovaSeq2.FCA/1277/{sample}_L00{l}_R1_001.fastq.gz",
        fq2 = "/staging/leuven/stg_00096/home/dkresa/ResolveOME/single_cell/230623.NovaSeq2.FCA/1277/{sample}_L00{l}_R2_001.fastq.gz",
        ref = reference,
        report_files = multiqc_output
    threads: 24
    output:
        bam = "raw_bams/{sample}_L00{l}.bam",
        task_done = "raw_bams/reports/bwa_mem2/{sample}_L00{l}.txt"
    shell:
        """
        /lustre1/project/stg_00096/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem \
        -t {threads} \
        -M \
        -R "@RG\\tID:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}" {input.ref} \
        {input.fq1} {input.fq2} \
        | samtools view -Sb -> {output.bam}

        touch {output.task_done}
        """

# Rule for sorting BAM files by queryname
rule sort_bam_files:
    input:
        bam = "raw_bams/{sample}_L00{l}.bam",
        check_input = bwa_mem2_report
    output:
        sorted_bam = "raw_bams/{sample}_L00{l}.sorted.bam",
        task_done = "raw_bams/reports/sort_bam_files/{sample}_L00{l}.txt"
    threads: 24
    shell:
        """
        java -jar $EBROOTPICARD/picard.jar SortSam \
        -I {input.bam} \
        -O {output.sorted_bam} \
        --SORT_ORDER queryname

        touch {output.task_done}
        """

# Rule for marking duplicates
rule mark_duplicates:
    input:
        bam = "raw_bams/{sample}_L00{l}.sorted.bam",
        check_input = sort_bam_files_report
    output:
        marked_bam = "raw_bams/{sample}_L00{l}.sorted.marked.bam",
        task_done = "raw_bams/reports/mark_duplicates/{sample}_L00{l}.txt"
    params:
        gatk = GATK_path,
        stage = stage_dir,
        lustre = lustre_dir
    threads: 24
    shell:
        """        
        singularity run --nv \
        -B {params.stage} \
        -B {params.lustre} \
        {params.gatk} gatk MarkDuplicatesSpark \
        -I {input.bam} \
        -O {output.marked_bam} \
        -M {output.marked_bam}.metrics.txt

        touch {output.task_done}
        """

# Rule for base quality score recalibration
rule base_recalibrator:
    input:
        marked_bam = "raw_bams/{sample}_L00{l}.sorted.marked.bam",
        ref = reference,
        dbsnp = known_sites,
        check_input = mark_duplicates_report
    output:
        data_table = "recalibrated/{sample}_L00{l}.recal_data.table",
        task_done = "recalibrated/reports/base_recalibrator/{sample}_L00{l}.txt"
    params:
        gatk = GATK_path,
        stage = stage_dir,
        lustre = lustre_dir
    threads: 2
    shell:
        """
        singularity run --nv \
        -B {params.stage} \
        -B {params.lustre} \
        {params.gatk} gatk BaseRecalibrator \
        -I {input.marked_bam} \
        -R {input.ref} \
        --known-sites {input.dbsnp} \
        -O {output.data_table}

        touch {output.task_done}
        """

# Rule for applying base quality score recalibration
rule apply_bqsr:
    input:
        marked_bam = "raw_bams/{sample}_L00{l}.sorted.marked.bam",
        recalibration_table = "recalibrated/{sample}_L00{l}.recal_data.table",
        ref = reference,
        check_input = base_recalibrator_report
    output:
        output_bam = "recalibrated/{sample}_L00{l}.sorted.marked.recal.bam",
        task_done = "recalibrated/reports/apply_bqsr/{sample}_L00{l}.txt"
    params:
        gatk = GATK_path,
        stage = stage_dir,
        lustre = lustre_dir
    threads: 2
    shell:
        """
        singularity run --nv \
        -B {params.stage} \
        -B {params.lustre} \
        {params.gatk} gatk ApplyBQSR \
        -R {input.ref} \
        -I {input.marked_bam} \
        --bqsr-recal-file {input.recalibration_table} \
        -O {output.output_bam}

        touch {output.task_done}
        """
# Rule for merging marked and recalibrated BAMs
rule merge_sorted_bam:
    input:
        bam_sorted_L1 = "recalibrated/{sample}_L001.sorted.marked.recal.bam",
        bam_sorted_L2 = "recalibrated/{sample}_L002.sorted.marked.recal.bam",
        bam_sorted_L3 = "recalibrated/{sample}_L003.sorted.marked.recal.bam",
        bam_sorted_L4 = "recalibrated/{sample}_L004.sorted.marked.recal.bam",
        report_L1 = "recalibrated/reports/apply_bqsr/{sample}_L001.txt",
        report_L2 = "recalibrated/reports/apply_bqsr/{sample}_L002.txt",
        report_L3 = "recalibrated/reports/apply_bqsr/{sample}_L003.txt",
        report_L4 = "recalibrated/reports/apply_bqsr/{sample}_L004.txt"
    output:
        merged_bam = "bam/{sample}.bam",
        task_done = "bam/reports/merge_sorted_bam/{sample}.txt"
    threads: 24
    shell:
        """
        sambamba merge -t {threads} {output.merged_bam} {input.bam_sorted_L1} {input.bam_sorted_L2} {input.bam_sorted_L3} {input.bam_sorted_L4}
        samtools index {output.merged_bam}

        touch {output.task_done}
        """

# Rule for haplotype calling
rule HaplotypeCaller:
    input:
        ref = reference,
        merged_bam = "bam/{sample}.bam",
        check_input = merge_sorted_bam_report,
        chrom = "contigs/chr{n}"
    output:
        gvcfs = "gvcf/{sample}/{sample}_chr{n}.g.vcf.gz",
        task_done = "reports/haplotype_calling/{sample}/{sample}_chr{n}.txt"
    params:
        gatk = GATK_path,
        stage = stage_dir,
        lustre = lustre_dir
    log:
       out = "logs/HaplotypeCaller/{sample}/{sample}_chr{n}.out",
       err = "logs/HaplotypeCaller/{sample}/{sample}_chr{n}.err"
    threads: 2
    shell:
        """
        singularity run --nv \
        -B {params.stage} \
        -B {params.lustre} \
        {params.gatk} gatk --java-options "-Xmx15g" HaplotypeCaller \
        -R {input.ref} \
        -I {input.merged_bam} \
        --intervals chr{wildcards.n} \
        -O {output.gvcfs} \
        -ERC GVCF \
        1> {log.out} 2> {log.err}

        touch {output.task_done}
        """

# Rule for merging GVCF files
rule GatherVcfs:
    input:
        ref = reference,
        check_input = haplotype_calling_report
    output:
        merged_gvcfs = "GatherVcfs/{sample}.g.vcf.gz",
        gvcf_index = "GatherVcfs/{sample}.g.vcf.gz.tbi",
        task_done = "reports/GatherVcfs/{sample}.txt"
    params:
        gatk = GATK_path,
        stage = stage_dir,
        lustre = lustre_dir,
        gvcfs = "{sample}_gvcf.list"
    threads: 2
    log:
       out = "logs/GatherVcfs/{sample}.out",
       err = "logs/GatherVcfs/{sample}.err"
    shell:
        """
        ls -1v gvcf/{wildcards.sample}/*.g.vcf.gz > {wildcards.sample}_gvcf.list

        singularity run --nv \
        -B {params.stage} \
        -B {params.lustre} \
        {params.gatk} gatk GatherVcfs \
        -I {params.gvcfs} \
        -O {output.merged_gvcfs} \
        -RI \
        --CREATE_INDEX true \
        1> {log.out} 2> {log.err}

        tabix -p vcf {output.merged_gvcfs}

        touch {output.task_done}
        """

# Rule for GVCFs combining (single-cell + bulk data)
rule GenomicsDBImport:
    input:
        ref = reference,
        check_input = GatherVcfs_report,
        chrom = "contigs/chr{n}"
    output:
        database_dir = directory("database/chr{n}_gdb"),
        task_done = "reports/GenomicsDBImport/chr{n}.txt"
    params:
        gatk = GATK_path,
        stage = stage_dir,
        lustre = lustre_dir,
        tmp = scratch_dir,
        sample_map = "sample.map"
    threads: 2
    log:
       out = "logs/GenomicsDBImport/chr{n}.out",
       err = "logs/GenomicsDBImport/chr{n}.err"
    shell:
        """
        for file in GatherVcfs/*.g.vcf.gz; do echo -e "$(basename "$file" .g.vcf.gz)\t$file"; done > sample.map

        singularity run --nv \
        -B {params.stage} \
        -B {params.lustre} \
        {params.gatk} gatk --java-options "-Xmx15g" GenomicsDBImport \
        --genomicsdb-workspace-path {output.database_dir} \
        -R {input.ref} \
        --sample-name-map {params.sample_map} \
        --intervals chr{wildcards.n} \
        --genomicsdb-shared-posixfs-optimizations true \
        1> {log.out} 2> {log.err}

        touch {output.task_done}
        """

# Rule for joint genotyping
rule GenotypeGVCFs:
    input:
        ref = reference,
        check_input = GenomicsDBImport_report,
        database = "database/chr{n}_gdb"
    output:
        vcf_chr = "vcf/chr{n}.vcf",
        task_done = "reports/GenotypeGVCFs/chr{n}.txt"
    params:
        gatk = GATK_path,
        stage = stage_dir,
        lustre = lustre_dir
    threads: 2
    log:
       out = "logs/GenotypeGVCFs/chr{n}.out",
       err = "logs/GenotypeGVCFs/chr{n}.err"
    shell:
        """
        singularity run --nv \
        -B {params.stage} \
        -B {params.lustre} \
        {params.gatk} gatk --java-options "-Xmx15g" GenotypeGVCFs \
        -R {input.ref} \
        --genomicsdb-shared-posixfs-optimizations true \
        -V gendb://{input.database} \
        --intervals chr{wildcards.n} \
        -O {output.vcf_chr} \
        1> {log.out} 2> {log.err}

        touch {output.task_done}
        """

# Rule for merging GVCF files
rule GatherVcfs2:
    input:
        ref = reference,
        check_input = GenotypeGVCFs_report
    output:
        merged_vcf = "vcf/merged.vcf",
        task_done = "reports/GatherVcfs2/merged.txt"
    params:
        gatk = GATK_path,
        stage = stage_dir,
        lustre = lustre_dir,
        vcfs = "merged_vcf.list"
    threads: 2
    log:
       out = "logs/GatherVcfs2/merged.out",
       err = "logs/GatherVcfs2/merged.err"
    shell:
        """
        ls -1v vcf/*.vcf > merged_vcf.list

        singularity run --nv \
        -B {params.stage} \
        -B {params.lustre} \
        {params.gatk} gatk GatherVcfs \
        -I {params.vcfs} \
        -O {output.merged_vcf} \
        -RI \
        --CREATE_INDEX true \
        1> {log.out} 2> {log.err}

        touch {output.task_done}
        """

# Rule for variant filtration
rule VariantFiltration:
    input:
        ref = reference,
        vcf = "vcf/merged.vcf",
        check_input = GatherVcfs2_report
    output:
        filtered_vcf = "vcf/merged_filtered.vcf",
        task_done = "reports/VariantFiltration/merged_filtered.txt"
    params:
        gatk = GATK_path,
        stage = stage_dir,
        lustre = lustre_dir
    threads: 2
    log:
       out = "logs/VariantFiltration/filtered.out",
       err = "logs/VariantFiltration/filtered.err"
    shell:
        """
        singularity run --nv \
        -B {params.stage} \
        -B {params.lustre} \
        {params.gatk} gatk --java-options "-Xmx15g" VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -O {output.filtered_vcf} \
        --filter-expression "QD < 2.0" \
        --filter-expression "MQ < 40.0" \
        --filter-expression "HaplotypeScore > 13.0" \
        --filter-expression "FS > 60.0" \
        --filter-expression "MQRankSum < -12.5" \
        --filter-expression "ReadPosRankSum < -8.0" \
        --filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
        --filter-expression "DP < 5" \
        --filter-expression "QUAL < 30.0" \
        --filter-expression "QUAL >= 30.0 && QUAL < 50.0" \
        --filter-expression "SOR > 4.0" \
        --filter-name "SNP_LowQualityDepth" \
        --filter-name "SNP_MappingQuality" \
        --filter-name "SNP_StrandBias" \
        --filter-name "SNP_HaplotypeScoreHigh" \
        --filter-name "SNP_MQRankSumLow" \
        --filter-name "SNP_ReadPosRankSumLow" \
        --filter-name "SNP_HardToValidate" \
        --filter-name "SNP_LowCoverage" \
        --filter-name "SNP_VeryLowQual" \
        --filter-name "SNP_LowQual" \
        --filter-name "SNP_SOR" \
        -cluster 3 \
        -window 10 \
        1> {log.out} 2> {log.err}

        touch {output.task_done}
        """

# Rule for quality control after mapping (CollectWgsMetrics)
rule CollectWgsMetrics:
    input:
        check_input = merge_sorted_bam_report,
        merged_bam_markdup = "bam/{sample}.bam",
        ref = reference
    output:
        wgs_metrics = "QC/2_postmap/{sample}_wgs_metrics.txt",
        task_done = "reports/CollectWgsMetrics/{sample}_wgs_metrics.txt"
    threads: 1
    log:
       out = "logs/CollectWgsMetrics/{sample}.out",
       err = "logs/CollectWgsMetrics/{sample}.err"
    shell:
        """
        java -jar $EBROOTPICARD/picard.jar CollectWgsMetrics \
        I={input.merged_bam_markdup} \
        O={output.wgs_metrics} \
        R={input.ref} \
        1> {log.out} 2> {log.err}

        touch {output.task_done}
        """

# Rule for quality control after mapping (CollectAlignmentSummaryMetrics)
rule CollectAlignmentSummaryMetrics:
    input:
        check_input = merge_sorted_bam_report,
        merged_bam_markdup = "bam/{sample}.bam",
        ref = reference
    output:
        alignment_summary_metrics = "QC/2_postmap/{sample}_alignment_summary_metrics.txt",
        task_done = "reports/CollectAlignmentSummaryMetrics/{sample}_alignment_summary_metrics.txt"
    threads: 1
    log:
       out = "logs/CollectAlignmentSummaryMetrics/{sample}.out",
       err = "logs/CollectAlignmentSummaryMetrics/{sample}.err"
    shell:
        """
        java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
        R={input.ref} \
        I={input.merged_bam_markdup} \
        O={output.alignment_summary_metrics} \
        1> {log.out} 2> {log.err}

        touch {output.task_done}
        """
