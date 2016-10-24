#!/usr/bin/perl -w
use strict;

my $dir = '../';
my $threads = 32;
my @organisms = ('yeast', 'mouse');
my @folds = ('25x');
my @strands = ('sample_01', 'sample_02');
my %ref_set = ('yeast', './yeast_cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.500.fa',
               'mouse', './mouse_cdna/Mus_musculus.GRCm38.cdna.all.500.fa');
my %index_set = ('yeast', './yeast_cdna/Saccharomyces_cerevisiae_cdna_500_bowtie2',
                 'mouse', './mouse_cdna/Mus_musculus_cdna_500_bowtie2');

chdir $dir;
foreach my $organism (@organisms){
    foreach my $strand (@strands) {
        foreach my $fold (@folds) {
            my $ref = $ref_set{$organism};
            my $index = $index_set{$organism};
            my $folder = "${organism}_${fold}";
            my $max_read_length = 100;
            system("/root/anaconda3/bin/python3 ./scripts/fasta_to_fastq.py ./$folder/${strand}_1.fasta > ./$folder/${strand}_1.fastq");
            system("/root/anaconda3/bin/python3 ./scripts/fasta_to_fastq.py ./$folder/${strand}_2.fasta > ./$folder/${strand}_2.fastq");
            my $read1 = "./$folder/${strand}_1.fastq";
            my $read2 = "./$folder/${strand}_2.fastq";

            system("mkdir -p ./$folder/bowtie2");
            my $mean = 250;
            my $sd = 25;

            my $min = $mean - 2 * $sd;
            my $max = $mean + 2 * $sd;
            
            my $start_time = time;
            system("bowtie2 --no-mix --no-discordant --end-to-end -a -p $threads -I $min -X $max --gbar $max_read_length -x $index -1 $read1 -2 $read2 > ./$folder/bowtie2/$strand.sam 2> ./$folder/bowtie2/${strand}_summary.txt");
            printf STDERR "- bowtie2 ${folder}_${strand}: %d sec\n", time - $start_time;

            $start_time = time;
            system("samtools view -@ $threads -b $folder/bowtie2/$strand.sam | samtools sort -@ $threads -n -o $folder/bowtie2/$strand.bam - > /dev/null 2> /dev/null");
            printf STDERR "- samtools view/sort ${folder}_${strand}: %d sec\n", time - $start_time;
            
            $start_time = time;
            system("express -o ./$folder/bowtie2_express_${strand} -m $mean -s $sd --fr-stranded $ref $folder/bowtie2/$strand.bam > ./log/bowtie2_express_${folder}_$strand.txt 2> /dev/null");
            printf STDERR "- express ${folder}_${strand}: %d sec\n", time - $start_time;
        }
    }
}
