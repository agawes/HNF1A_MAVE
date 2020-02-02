# HNF1A_MAVE

Preliminary analysis of HNF1A MAVE sequencing data

1. "Pipeline" - concatenate appropriate FASTQ files from different lanes, trim sequencing adapters, map with bwa mem, samtools fixmates, and run custom script for counting HNF1A codons (tri-nucleotides): [script](script.sh)

2. Custom script to count HNF1A tri-nucleotides, and compare them to the MAVE design file: [script](script.sh)

3. Basic QC plots - how many design and non-design codons we see at different depth of coverage: [script](script.sh)
