#!/usr/bin/perl -w
use strict;

my $dir = '../';
my $threads = 32;
my @organisms = ('yeast', 'mouse');
my @folds = ('10x');
my @strands = ('sample_01', 'sample_02');
my %ref_set = ('yeast', './yeast_cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.500.fa',
               'mouse', './mouse_cdna/Mus_musculus.GRCm38.cdna.all.500.fa');
my %index_set = ('yeast', './yeast_cdna/Saccharomyces_cerevisiae_cdna_500_rsem_bowtie2',
                 'mouse', './mouse_cdna/Mus_musculus_cdna_500_rsem_bowtie2');

chdir $dir;
foreach my $organism (@organisms){
    foreach my $strand (@strands) {
        foreach my $fold (@folds) {
            my $index = $index_set{$organism};
            my $folder = "${organism}_${fold}";
            my $max_read_length = 100;
            
            #system("/root/anaconda3/bin/python3 ./scripts/fasta_to_fastq.py ./$folder/${strand}_1.fasta > ./$folder/${strand}_1.fastq");
            #system("/root/anaconda3/bin/python3 ./scripts/fasta_to_fastq.py ./$folder/${strand}_2.fasta > ./$folder/${strand}_2.fastq");
            
            my $read1 = "./$folder/${strand}_1.fastq";
            my $read2 = "./$folder/${strand}_2.fastq";

            my $mean = 250;
            my $sd = 25;
            
            system("mkdir -p ./$folder/bowtie2_rsem_$strand");

            my $start_time = time;
            system("rsem-calculate-expression --forward-prob 1 --fragment-length-mean $mean --fragment-length-sd $sd -p $threads --bowtie2 --paired-end --estimate-rspd --time $read1 $read2 $index ./$folder/bowtie2_rsem_$strand/bowtie2_rsem_$strand > ./log/rsem_${folder}_$strand.txt 2> ./log/rsem_${folder}_$strand.err");
            printf STDERR "- rsem ${folder}_${strand}: %d sec\n", time - $start_time;
=comment
            system("samtools sort -@ $threads -n -o ./$folder/bowtie2_rsem_$strand/bowtie2_rsem_$strand.transcript.sorted.bam ./$folder/bowtie2_rsem_$strand/bowtie2_rsem_$strand.transcript.bam ");
            system("express -o ./$folder/bowtie2_rsem_$strand/express_$strand -m $mean -s $sd --fr-stranded $ref_set{$organism} ./$folder/bowtie2_rsem_$strand/bowtie2_rsem_$strand.transcript.sorted.bam > ./log/rsem_param_bowtie2_express_${folder}_$strand.out 2> ./log/rsem_param_bowtie2_express_${folder}_$strand.err");
=cut
        }
    }
}
