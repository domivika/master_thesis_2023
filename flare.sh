# we start by combining all of fastq reads
cat fastq/*.fastq > combined_samples.fastq

flair align -t 24 \
-g "/path/to/Homo_sapiens_assembly38_masked.fasta" \
-r "/path/to/combined_samples.fastq" \
-o "combined_samples.flair.aligned"

flair correct -t 24 \
-g "/path/to/bwa-mem2/Homo_sapiens_assembly38_masked.fasta" \
--gtf "/path/to/gencode.v46.basic.annotation.gtf" \
-q combined_samples.flair.aligned.bed \
-o combined_samples.flair

flair collapse -t 24 \
-g "/path/to/Homo_sapiens_assembly38_masked.fasta" \
--gtf "/path/to/gencode.v46.basic.annotation.gtf" \
-q combined_samples.flair_all_corrected.bed \
-r "/path/to/combined_samples.fastq" \
--annotation_reliant generate --generate_map \
--check_splice --stringent \
--output combined_samples.flair.collapse

flair quantify -t 24 \
-r "/path/to/sample_manifest.tsv" \
-i combined_samples.flair.collapse.isoforms.fa \
--generate_map --isoform_bed combined_samples.flair.collapse.isoforms.bed \
--stringent --check_splice

     