setwd('/Users/Hsieh/Documents/Projects/Transcript Quantification/')

al_nr = read.table('mouse_data/answer_label_mouse_trinity_10x.tsv', header=F, sep='\t', stringsAsFactors=F)
al_split = read.table('mouse_data/answer_label_mouse_trinity_10x_answer_label_split.tsv', header=F, sep='\t', stringsAsFactors=F)

al_split_target = al_split[grep('split', al_split$V2), ]
rownames(al_nr) = al_nr$V1
rownames(al_split_target) = al_split_target$V1
al_nr_target = subset(al_nr, al_nr$V1 %in% al_split_target$V1)

al_nr_target = al_nr_target[sort(al_nr_target$V1), ]
al_split_target = al_split_target[sort(al_split_target$V1),]
 
answer = read.table('mouse_10x/mouse_10x_count_answer.tsv', header=T, sep='\t')
answer_target = subset(answer, answer$target_id %in% al_nr_target$V1)
answer_target = answer_target[order(answer_target$target_id, decreasing=F), ]
rownames(answer_target) = answer_target$target_id

nr_count = read.table('mouse_10x/bowtie2_rsem_sample_01_trinity_answer_label/bowtie2_rsem.isoforms.results', header=T, sep='\t')
rownames(nr_count) = nr_count$transcript_id
al_nr_target$expected_count = nr_count[al_nr_target$V2, 'expected_count']
split_count = read.table('mouse_10x/bowtie2_rsem_sample_01_trinity_answer_label_split_nr/bowtie2_rsem.isoforms.results', header=T, sep='\t')
rownames(split_count) = split_count$transcript_id
al_split_target$expected_count = split_count[al_split_target$V2, 'expected_count']
al_split_target$answer = answer_target[al_split_target$V1, 'sample_01']

reference = read.table('mouse_10x/bowtie2_rsem_sample_01/bowtie2_rsem_sample_01.isoforms.results', header=T, sep='\t')
reference_target = subset(reference, reference$transcript_id %in% answer_target$target_id)
reference_target = reference_target[order(reference_target$transcript_id, decreasing=F), ]
al_split_target$reference = reference_target$expected_count
