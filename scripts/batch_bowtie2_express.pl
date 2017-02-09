#!/usr/bin/perl -w
use strict;

my $dir = '/home/dn070017/projects/transcript_quatification';
my $threads = 32;
#my @organisms = ('yeast', 'mouse');
my @organims = ('mouse');
#my @targets = ('reference', 'trinity');
my @targets = ('reference');
#my @folds = ('10x', '25x');
my @folds = ('10x');
my @strands = ('sample_01', 'sample_02');

my %ref_fasta = ('yeast_10x_reference', './yeast_data/Saccharomyces_cerevisiae.R64-1-1.cdna.all.500.fa',
                 'mouse_10x_reference', './mouse_data/Mus_musculus.GRCm38.cdna.all.500.fa',
                 'yeast_10x_trinity', './yeast_data/yeast_trinity_10x.fasta',
                 'mouse_10x_trinity', './mouse_data/mouse_trinity_10x.fasta');

my %ref_index = ('yeast_10x_reference', './yeast_data/Saccharomyces_cerevisiae_cdna_500_bowtie2',
                 'mouse_10x_reference', './mouse_data/Mus_musculus_cdna_500_bowtie2',
                 'yeast_10x_trinity', './yeast_data/yeast_trinity_10x_bowtie2',
                 'mouse_10x_trinity', './mouse_data/mouse_trinity_10x_bowtie2');

chdir $dir;

foreach my $organism (@organisms) {
    foreach my $target (@targets) {
        foreach my $strand (@strands) {
            foreach my $fold (@folds) {
                my $folder = "${organism}_${fold}";
                my $suffix = "${organism}_${fold}_${target}";
                my $fasta = $ref_fasta{$suffix};
                my $index = $ref_index{$suffix};
                
                open TIME, ">./time_log/bowtie2_express_${folder}_${strand}_$target.txt" or die 
                           "can not open ./time_log/bowtie2_express_${folder}_${strand}_$target.txt";
            
                my $read1 = "./$folder/${strand}_1.fastq";
                my $read2 = "./$folder/${strand}_2.fastq";

                system("mkdir -p ./$folder/bowtie2_$target");
            
                my $mean = 250;
                my $sd = 25;

                my $min = $mean - 2 * $sd;
                my $max = $mean + 2 * $sd;
                my $max_read_length = 100;
                
                my $start_time = time;
                system("bowtie2 --no-mix --no-discordant --end-to-end -a -p $threads -I $min -X $max --gbar $max_read_length -x $index -1 $read1 -2 $read2 > ./$folder/bowtie2_$target/$strand.sam 2> ./$folder/bowtie2_$target/${strand}_summary.txt");
                printf TIME "- bowtie2 ${folder}_${strand}_$target: %d sec\n", time - $start_time;

                $start_time = time;
                system("samtools view -@ $threads -b $folder/bowtie2_$target/$strand.sam | samtools sort -@ $threads -n -o $folder/bowtie2_$target/$strand.bam - > /dev/null 2> /dev/null");
                printf TIME "- samtools view/sort ${folder}_${strand}_$target: %d sec\n", time - $start_time;
            
                $start_time = time;
                system("mkdir -p ./$folder/bowtie2_express_${strand}_$target");
                system("express -o ./$folder/bowtie2_express_${strand}_$target -m $mean -s $sd --fr-stranded $fasta $folder/bowtie2_$target/$strand.bam > ./$folder/bowtie2_express_${strand}_$target/express.out 2> /dev/null");
                printf TIME "- express ${folder}_${strand}_$target: %d sec\n", time - $start_time;
            }
        }
    }
}
