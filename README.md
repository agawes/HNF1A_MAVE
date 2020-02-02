# HNF1A_MAVE

Preliminary analysis of HNF1A MAVE sequencing data

1. "Pipeline" - concatenate appropriate FASTQ files from different lanes, trim sequencing adapters, map with bwa mem, samtools fixmates, and run custom script for counting HNF1A codons (tri-nucleotides): [script](hnf1a_mave_2020.pipeline.sh)

2. Custom script to count HNF1A tri-nucleotides, and compare them to the MAVE design file: [script](get_HNF1A_codon_summary.pl)

3. Basic QC plots - how many design and non-design codons we see at different depth of coverage: [script](QC_plots.R)
