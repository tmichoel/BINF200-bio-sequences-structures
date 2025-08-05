#!/usr/bin/bash

# Set the data folder
sampledir=../data/PacBio_VSG
cd $sampledir

# Database file
dbf=TREU927-v26_VSGTranscripts/TREU927-v26_VSGTranscripts.fasta

# Run for all sample files
for qf in filtered_reads/PacBio_VSG_filtered_reads_*.fasta
do 
    # Sample name
    sample=`echo "$qf" | cut -d'_' --fields=6-8 | cut -d'.' -f1`
    echo $sample
    # Blast output file name
    of=blastn/PacBio_VSG_filtered_reads_blastn_$sample.txt
    # Blast command
    blastn  -query $qf -db $dbf -out $of  -max_target_seqs 1 -outfmt '6 qseqid sseqid score bitscore evalue qlen slen length sstart send nident mismatch gaps positive'
    # ORF file name
    orf=orf/PacBio_VSG_filtered_reads_ORF_$sample.fasta
    # GetORF command
    /usr/bin/getorf -sequence $qf -outseq $orf -minsize 1200 -find 3 -reverse N
done

