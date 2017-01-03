library(PerformanceAnalytics)

setwd('/Users/Hsieh/Documents/Projects/Transcript_Quantification/')

organism = 'yeast'
#organism = 'mouse'

fold = '10x'

target = paste0(organism, '_', fold)
answer = read.table(paste0(target, '/', target, '_count_answer.tsv'), header=T, sep='\t', 
                    stringsAsFactors=F)
answer = answer[order(answer$target_id), ]
colnames(answer) = c('target_id', 'answer_1', 'answer_2')

#path_1 = paste0(target, '/rsem_param_bowtie2_express_sample_01/results.xprs')
#path_2 = paste0(target, '/rsem_param_bowtie2_express_sample_02/results.xprs')
#path_1 = paste0(target, '/bowtie2_express_sample_01/results.xprs')
#path_2 = paste0(target, '/bowtie2_express_sample_02/results.xprs')
path_1 = paste0(target, '/bowtie2_express_sample_01_trinity/results.xprs')
path_2 = paste0(target, '/bowtie2_express_sample_02_trinity/results.xprs')
raw_data_1 = read.table(path_1, header=T, sep='\t', stringsAsFactors=F)
raw_data_1 = raw_data_1[order(raw_data_1$target_id), c('target_id', 'eff_counts')]
raw_data_2 = read.table(path_2, header=T, sep='\t', stringsAsFactors=F)
raw_data_2 = raw_data_2[order(raw_data_2$target_id), c('target_id', 'eff_counts')]

path_3 = paste0(organism, '_cdna/', organism, '_trinity_', fold, '_completness.tsv')
completness = read.table(path_3, header=F, stringsAsFactors=F, sep='\t') 
colnames(completness) = c('target_id', 'set_id', 'orientation', 'target_length', 'set_length',
                          'target_alignment', 'set_alignment', 'target_coverage', 'set_coverage')
completness = completness[order(completness$target_id),]
answer = subset(answer, answer$target_id %in% completness$target_id)
rownames(raw_data_1) = raw_data_1$target_id
rownames(raw_data_2) = raw_data_2$target_id

merge_all = cbind(answer[,c(2,3)], raw_data_1[completness$set_id, 'eff_counts'], 
                  raw_data_2[completness$set_id, 'eff_counts'])
colnames(merge_all) = c('answer_1', 'answer_2', 'express_1', 'express_2')
rownames(merge_all) = answer$target_id

merge_all_full = subset(merge_all, completness$target_coverage >= 0.9 & completness$set_coverage >= 0.9)

#png(filename=paste0('figures/rsem_param_express_', organism, '_', fold, '_total.png'), width=800, height=750)    
png(filename=paste0('figures/express_', organism, '_', fold, '_total.png'), width=800, height=750)    
chart.Correlation(merge_all_full, histogram=TRUE, method="pearson", pch=20)
dev.off()

blast_summary = read.table(paste0(organism, '_cdna/blastn_pairwise_', organism, 
                           '_cdna_100_ambiguous.tsv'), header=T, sep='\t',
                           stringsAsFactors=F)
unique_set = subset(blast_summary, blast_summary$ambiguous_sites == 0)
ambiguous_set = subset(blast_summary, blast_summary$ambiguous_sites > 0)

# print blast summary
print(paste0('total: ', length(blast_summary$target_id)))
print(paste0('unique: ', length(unique_set$target_id)))
print(paste0('ambiguous: ', length(ambiguous_set$target_id)))

# plot result
merge_unique = subset(merge_all, !rownames(merge_all) %in% ambiguous_set$target_id)
#png(filename=paste0('figures/rsem_param_express_', organism, '_', fold, '_unique.png'), width=800, height=750)
png(filename=paste0('figures/express_', organism, '_', fold, '_unique.png'), width=800, height=750)
chart.Correlation(merge_unique, histogram=TRUE, method="pearson", pch=20)
dev.off()

merge_ambiguous = subset(merge_all, rownames(merge_all) %in% ambiguous_set$target_id)
#png(filename=paste0('figures/rsem_param_express_', organism, '_', fold, '_ambiguous.png'), width=800, height=750)
png(filename=paste0('figures/express_', organism, '_', fold, '_ambiguous.png'), width=800, height=750)
chart.Correlation(merge_ambiguous, histogram=TRUE, method="pearson", pch=20)
dev.off()



