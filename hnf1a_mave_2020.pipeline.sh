cd /well/mccarthy/production/genetic-screens/data/HNF1A_MAVE/TWIST_Jan2020/bsg-ftp.well.ox.ac.uk/200103_NB502094_0176_AHVHL5BGXC/

perl /well/mccarthy/production/rna-seq/code/pipeline_perl/0.get_RG_info.pl .

## 20 samples - not sure what is what...
## each was sequenced across 4 seq lanes

## moved raw fastq's to raw_fastq folder, and RG_info to top folder; rest can be deleted

## concat appropriate fastq's using RG_info
perl concat_fastq.pl &

## trim adapters
mkdir trim_fastq
for I in concat_fastq/*.1.fastq.gz; do
 base=`basename $I .1.fastq.gz`
  printf "$base\n"
  echo "/well/got2d/agata/bin/trim_galore/trim_galore -o trim_fastq/ --paired $I concat_fastq/$base.2.fastq.gz" > trim_fastq/$base.trim.sh
  qsub -P mccarthy.prjc -q short.qc -V -cwd -N trim_$base -e trim_fastq/$base.trim.err -o trim_fastq/$base.trim.out trim_fastq/$base.trim.sh
  printf "$base done\n"
done

##### bwa mem
mkdir bwa_mem
for I in trim_fastq/*.1_val_1.fq.gz; do
        base=`basename $I .1_val_1.fq.gz`
        printf "$I\t$base\n"
        echo "/apps/well/bwa/0.7.12/bwa mem ../MAVE_SALT_1/reference/HNF1A_whole_construct.fna $I trim_fastq/$base.2_val_2.fq.gz > bwa_mem/$base.bwa_mem.sam" > bwa_mem/$base.bwa_mem.sh
        qsub -P mccarthy.prjc -q short.qc -V -cwd -N bwa_$base -e bwa_mem/$base.mem.err -o bwa_mem/$base.mem.out bwa_mem/$base.bwa_mem.sh
done

### run samtools fixmates - fill in mate coordinates, ISIZE and mate related flags from a name-sorted alignments
mkdir fixmate_bam
for I in bwa_mem/*.bwa_mem.sam; do
        base=`basename $I .bwa_mem.sam`
        printf "$I\t$base\n"
        echo "samtools fixmate -r -O sam $I fixmate_bam/$base.bwa_mem.fixmate.sam" > fixmate_bam/$base.fixmate.sh
        qsub -P mccarthy.prjc -q short.qc -V -cwd -N fixmate_$base -e fixmate_bam/$base.fixmate.err -o fixmate_bam/$base.fixmate.out fixmate_bam/$base.fixmate.sh
done

mkdir codon_summary
for I in fixmate_bam/*.bwa_mem.fixmate.sam; do
        base=`basename $I .bwa_mem.fixmate.sam`
        printf "$I\t$base\n"
        echo "perl /well/mccarthy/production/genetic-screens/code/get_HNF1A_codon_summary.pl -s $I > codon_summary/$base.summary" > codon_summary/$base.codons.sh
        qsub -P mccarthy.prjc -pe shmem 12 -q short.qc -V -cwd -N codons_$base -e codon_summary/$base.codons.err -o codon_summary/$base.codons.out codon_summary/$base.codons.sh
done