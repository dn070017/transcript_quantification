#!/usr/bin/perl -w
use strict;

my $dir = '../';
my $threads = 32;
my @organisms = ('yeast', 'mouse');
my @folds = ('10x');
my @kmers = (1, 2, 3);
my @strands = ('sample_01', 'sample_02');

chdir $dir;
foreach my $organism (@organisms){
    foreach my $kmer (@kmers) {
        foreach my $fold (@folds) {
            my $folder = "${organism}_${fold}";
            
            #system("/root/anaconda3/bin/python3 ./scripts/fasta_to_fastq.py ./$folder/${strand}_1.fasta > ./$folder/${strand}_1.fastq");
            #system("/root/anaconda3/bin/python3 ./scripts/fasta_to_fastq.py ./$folder/${strand}_2.fasta > ./$folder/${strand}_2.fastq");
            
            my $s1_read1 = "./$folder/$strands[0]_1.fastq";
            my $s1_read2 = "./$folder/$strands[0]_2.fastq";
            my $s2_read1 = "./$folder/$strands[1]_1.fastq";
            my $s2_read2 = "./$folder/$strands[1]_2.fastq";

            system("mkdir -p ./$folder/trinity_$kmer");
        
            my $mean = 250;
            my $sd = 25;
            my $max = $mean + 2 * $sd;

            my $start_time = time;
            system("Trinity --CPU $threads --max_memory 300G --seqType fq --left $s1_read1,$s2_read1 --right $s1_read2,$s2_read2 --group_pairs_distance $max --min_kmer_cov $kmer --SS_lib_type FR --output ./$folder/trinity_$kmer > ./log/trinity_${folder}_$kmer.out 2> ./log/trinity_${folder}_$kmer.err");
            printf STDERR "- trinity $folder: %d sec\n", time - $start_time;
            
            system("mv ./$folder/trinity_$kmer/Trinity.fasta ./${organism}_cdna/${folder}_$kmer.fasta");
        }
    }
}
