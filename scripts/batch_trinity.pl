#!/usr/bin/perl -w
use strict;

my $dir = '../';
my $threads = 32;
my @organisms = ('yeast', 'mouse');
my @folds = ('25x');

chdir $dir;

foreach my $organism (@organisms){
    foreach my $fold (@folds) {
        my $folder = "${organism}_${fold}";
        
        open TIME, ">./time_log/trinity_${folder}.txt" or die "can not open ./time_log/trinity_${folder}.txt";
            
        #system("/root/anaconda3/bin/python3 ./scripts/fasta_to_fastq.py ./$folder/${strand}_1.fasta > ./$folder/${strand}_1.fastq");
        #system("/root/anaconda3/bin/python3 ./scripts/fasta_to_fastq.py ./$folder/${strand}_2.fasta > ./$folder/${strand}_2.fastq");
            
        my $s1_read1 = "./$folder/sample_01_1.fastq";
        my $s1_read2 = "./$folder/sample_01_2.fastq";
        my $s2_read1 = "./$folder/sample_02_1.fastq";
        my $s2_read2 = "./$folder/sample_02_2.fastq";

        system("mkdir -p ./$folder/trinity");
        
        my $mean = 250;
        my $sd = 25;
        my $max = $mean + 2 * $sd;

        my $start_time = time;
        system("mkdir -p ./$folder/trinity");
        system("Trinity --CPU $threads --max_memory 300G --seqType fq --left $s1_read1,$s2_read1 --right $s1_read2,$s2_read2 --group_pairs_distance $max --SS_lib_type FR --output ./$folder/trinity > ./$folder/trinity/trinity.out 2> ./$folder/trinity/trinity.err");
        printf TIME "- trinity $folder: %d sec\n", time - $start_time;

        system("mv ./$folder/trinity/Trinity.fasta ./${organism}_data/${organism}_trinity_$fold.fasta");

        close TIME;
    }
}
