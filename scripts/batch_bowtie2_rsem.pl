#!/usr/bin/perl -w
use strict;

my $dir = '/home/dn070017/projects/transcript_quatification';
my $threads = 32;
my @organisms = ('yeast', 'mouse');
@organisms = ('mouse');
#my @targets = ('reference', 'trinity');
my @targets = ('trinity_answer_label_split_nr');
#my @folds = ('10x', '25x');
my @folds = ('10x');
my @strands = ('sample_01', 'sample_02');

my %ref_index = ('yeast_10x_reference', './yeast_data/Saccharomyces_cerevisiae_cdna_500_rsem_bowtie2',
                 'mouse_10x_reference', './mouse_data/Mus_musculus_cdna_500_rsem_bowtie2',
                 'yeast_10x_trinity', './yeast_data/yeast_trinity_10x_rsem_bowtie2',
                 'mouse_10x_trinity', './mouse_data/mouse_trinity_10x_rsem_bowtie2');

%ref_index = ('yeast_10x_reference', './yeast_data/Saccharomyces_cerevisiae_cdna_500_rsem_bowtie2',
              'mouse_10x_reference', './mouse_data/Mus_musculus_cdna_500_rsem_bowtie2',
              'yeast_10x_trinity_answer_label', './yeast_data/yeast_trinity_10x_answer_label_rsem_bowtie2',
              'mouse_10x_trinity_answer_label', './mouse_data/mouse_trinity_10x_answer_label_rsem_bowtie2',
              'mouse_10x_trinity_answer_label_split_nr', './mouse_data/mouse_trinity_10x_answer_label_split_nr_rsem_bowtie2');
chdir $dir;

foreach my $organism (@organisms) {
    foreach my $target (@targets) {
        foreach my $strand (@strands) {
            foreach my $fold (@folds) {
                my $folder = "${organism}_${fold}";
                my $suffix = "${organism}_${fold}_${target}";
                my $index = $ref_index{$suffix};
                
                open TIME, ">./time_log/bowtie2_rsem_${folder}_${strand}_$target.txt" or die 
                           "can not open ./time_log/bowtie2_rsem_${folder}_${strand}_$target.txt";
                
                my $read1 = "./$folder/${strand}_1.fastq";
                my $read2 = "./$folder/${strand}_2.fastq";

                my $mean = 250;
                my $sd = 25;
            
                system("mkdir -p ./$folder/bowtie2_rsem_${strand}_$target");

                my $start_time = time;
                system("rsem-calculate-expression --forward-prob 1 --fragment-length-mean $mean --fragment-length-sd $sd -p $threads --bowtie2 --paired-end --estimate-rspd --time $read1 $read2 $index ./$folder/bowtie2_rsem_${strand}_$target/bowtie2_rsem > ./$folder/bowtie2_rsem_${strand}_$target/rsem.out 2> ./$folder/bowtie2_rsem_${strand}_$target/rsem.err");
                printf TIME "- rsem ${folder}_${strand}: %d sec\n", time - $start_time;
            }
        }
    }
}
