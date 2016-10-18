#!/usr/bin/Rscript

suppressMessages(library('Biostrings'))
suppressMessages(library('polyester'))
suppressMessages(library('tools'))

args = commandArgs(trailingOnly=T)

# read command line arguments
if(length(args) != 4) {
    error_message = paste0('usage: simulation.R [target.fasta] ',
    '[percentage of diffential genes] [sequencing depth] ', 
    '[output prefix]\n')
    stop(error_message, call.=F)
} else {
    fasta_file = file_path_as_absolute(args[1])
    fasta_count = count_transcripts(fasta_file) 
    deg_percentage = as.numeric(args[2])
    depth = as.numeric(args[3])
    output_prefix = args[4]
}

fasta = readDNAStringSet(fasta_file)

# setting fold changes
cat('do you want to use user specific fold change setting? (yes/no) ')
load_fold_change = readLines('stdin', n=1)
fold_change_file = ''
if(load_fold_change == 'yes' || load_fold_change == 'y'){
    cat('please specify path for user specific fold change setting: ')
    fold_change_file = file_path_as_absolute(readLines('stdin', n=1))
    fold_change = read.table(fold_change_file, header=T, sep='\t')
} else {
    fold_change = data.frame(value=matrix(sample(seq(from=0.01, to=10.00, 
    by=0.01), fasta_count, replace=T), fasta_count))
    for(i in 1:fasta_count){
        if(runif(1, 0, 1) > deg_percentage) {
            fold_change[i, 1] = sample(seq(from=0.8, to=1.2, by=0.01), 1)
        }
    }
}

# print arguments and fold change setting
dir.create(output_prefix, recursive=T)
setwd(output_prefix)

cat(paste0('target.fasta: ', fasta_file, '\nfasta count: ', 
fasta_count, '\npercentage of diffential genes: ', deg_percentage, 
'\nsequncing depth: ', depth, '\noutput prefix: ', output_prefix,
'\nfold change setting: ', ifelse(fold_change_file != '', 
fold_change_file, 'randomly generated')), file='arguments.txt')

write.table(fold_change, 'fold_change.txt', quote=F, col.names=T, 
row.names=F)

# processing simulation
split = fasta_count / 250
if(split != fasta_count %/% 250) {
    split = fasta_count %/% 250 + 1
}

for(x in c(1:split)) {

    start = 250 * (x - 1) + 1
    end = 250 * x

    if(end > fasta_count) { end = fasta_count }

    cat(paste0('- processing sequences ', start, '-', end, '\n'))

    folder = paste('temporary_sequences/sequences', start, 'to', end, sep='_')
    subset_file = paste0(folder, '/', start, '_', end, ".fasta")

    dir.create(folder, recursive=T)
    
    subset_fasta = fasta[start:end]
    subset_fold_change = fold_change[start:end, 1]
    
    writeXStringSet(subset_fasta, subset_file)

    # ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
    # here all transcripts will have ~equal FPKM
    readspertx = round(depth * (width(subset_fasta) / 200))

    # simulation call
    # - paired end
    # - mean of fragment length: 250 
    # - sd of fragment length: 25
    # - error model: illumina5 (TruSeq), illumina4 (HiSeq SBS, HiSeq Cluster)
    # - strand specific: fr-stranded
    # - bias: rnaf (little different on 3' and 5' site specific bias) 
    simulate_experiment(subset_file, reads_per_transcript=readspertx, 
    num_reps=c(1,1), fold_changes=subset_fold_change, gc_bias=1, 
    error_model='illumina5', bias='rnaf', outdir=folder,
    strand_specific=T) 
}

system('cat */*/sample_01_1.fasta > sample_01_1.fasta')
system('cat */*/sample_01_2.fasta > sample_01_2.fasta')
system('cat */*/sample_02_1.fasta > sample_02_1.fasta')
system('cat */*/sample_02_2.fasta > sample_02_2.fasta')
