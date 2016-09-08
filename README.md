# Plasmodium falciparum transcription factor ChIP-seq

Data analysis of various ChIP-seq experiments.

Workflow is split up into three subworkflows:

1. Raw read processing, trimming, and alignment
2. Alignment processing and peak calling
3. Motif analysis and peak annotation

## Test workflow

A test workflow was designed to run on sampled data from `ap2i_{sample}_r1.fastq.gz` where `{sample}` is either chip, input, or igg:

```
seqtk sample -s 1213 ap2i_{sample}_r1.fastq.gz 500000 > sample_{sample}.fastq && gzip sample_{sample}.fastq
```

I sampled 500,000 reads from each library.
