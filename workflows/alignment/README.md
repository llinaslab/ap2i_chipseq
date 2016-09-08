# Alignment

Within this workflow we perform quality control measures on the raw read files, trim them of potential adapter contamination, align them to a reference genome, and filter that alignment to keep only high quality, non-duplicate, primary alignments.

## Workflow

1. Run fastqc
2. Trim reads using trimmomatic
3. Align reads using bwa mem & filter alignments using samtools
