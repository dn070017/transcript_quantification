setwd('/Users/Hsieh/Documents/Projects/Transcript Quantification/')

#organism = 'yeast'
organism = 'mouse'

fold = '10x'

prefix = paste0(organism, '_', fold)
answer = read.table(paste0(prefix, '/', prefix, '_count_answer.tsv'), header=T, sep='\t', 
                    stringsAsFactors=F)
answer = answer[order(answer$target_id), ]
rownames(answer) = answer$target_id
colnames(answer) = c('target_id', 'answer_1', 'answer_2')
meta_data$ans_1 = answer[meta_data$ref_id, 'answer_1']
meta_data$ans_2 = answer[meta_data$ref_id, 'answer_2']

# reference set
path_1 = paste0(prefix, '/bowtie2_rsem_sample_01/bowtie2_rsem_sample_01.isoforms.results')
path_2 = paste0(prefix, '/bowtie2_rsem_sample_02/bowtie2_rsem_sample_02.isoforms.results')
raw_data_1 = read.table(path_1, header=T, sep='\t', stringsAsFactors=F)
raw_data_1 = raw_data_1[order(raw_data_1$transcript_id), c('transcript_id', 'expected_count')]
raw_data_2 = read.table(path_2, header=T, sep='\t', stringsAsFactors=F)
raw_data_2 = raw_data_2[order(raw_data_2$transcript_id), c('transcript_id', 'expected_count')]
rownames(raw_data_1) = raw_data_1$transcript_id
rownames(raw_data_2) = raw_data_2$transcript_id
meta_data$ref_count_1 = raw_data_1[meta_data$ref_id, 'expected_count']
meta_data$ref_count_2 = raw_data_2[meta_data$ref_id, 'expected_count']
    
# trinity set
path_1 = paste0(prefix, '/bowtie2_rsem_sample_01_trinity_answer_label/bowtie2_rsem.isoforms.results')
path_2 = paste0(prefix, '/bowtie2_rsem_sample_02_trinity_answer_label/bowtie2_rsem.isoforms.results')
raw_data_1 = read.table(path_1, header=T, sep='\t', stringsAsFactors=F)
raw_data_1 = raw_data_1[order(raw_data_1$transcript_id), c('transcript_id', 'expected_count')]
raw_data_2 = read.table(path_2, header=T, sep='\t', stringsAsFactors=F)
raw_data_2 = raw_data_2[order(raw_data_2$transcript_id), c('transcript_id', 'expected_count')]
rownames(raw_data_1) = raw_data_1$transcript_id
rownames(raw_data_2) = raw_data_2$transcript_id
meta_data$set_count_1 = raw_data_1[meta_data$set_id, 'expected_count']
meta_data$set_count_2 = raw_data_2[meta_data$set_id, 'expected_count']

# assembly completness
for(target_category in unique(meta_data$category)) {
    target_meta = subset(meta_data, meta_data$category == target_category)
    figure_1 = ggplot(target_meta, aes(ans_1, set_count_1, shape='a', colour="#FF9999"))
    p_cor = round(cor(target_meta$ans_1, target_meta$set_count_1), 2)
    figure_1 = figure_1 + geom_point() + 
             ggtitle(paste0(organism, ' - ', target_category, ' (', p_cor, ')')) +
             scale_shape_discrete(solid=F, guide=F) + xlab('answer') + ylab(target) +
             scale_colour_discrete(guide=F) + xlim(0, quantile(target_meta$ans_1, 0.95)) + 
             ylim(0, quantile(target_meta$set_count_1, 0.95))
    
    target_meta$ans_1 = target_meta$ans_1 / target_meta$ref_length
    target_meta$set_count_1 = target_meta$set_count_1 / target_meta$set_length
    p_cor = round(cor(target_meta$ans_1, target_meta$set_count_1), 2)
    figure_2 = ggplot(target_meta, aes(ans_1, set_count_1, shape='a', colour="#FF9999"))
    figure_2 = figure_2 + geom_point() + 
             ggtitle(paste0('(normalized) ', organism, ' - ', target_category, ' (', p_cor, ')')) +
             scale_shape_discrete(solid=F, guide=F) + xlab('answer') + ylab(target) +
             scale_colour_discrete(guide=F) + xlim(0, quantile(target_meta$ans_1, 0.95)) + 
             ylim(0, quantile(target_meta$set_count_1, 0.95))
    
    png(paste0('figures/', organism, '_xy_plot_answer_label_', target_category, '_s1.png'), width=1000, height=400)
    multiplot(figure_1, figure_2, cols=2)
    dev.off()
}

# ambiguous/unique sites
for(target_category in unique(meta_data$ambiguous_type)) {
    target_meta = subset(meta_data, meta_data$ambiguous_type == target_category)
    figure_1 = ggplot(target_meta, aes(ans_1, set_count_1, shape='a', colour="#FF9999"))
    p_cor = round(cor(target_meta$ans_1, target_meta$set_count_1), 2)
    figure_1 = figure_1 + geom_point() + 
        ggtitle(paste0(organism, ' - ', target_category, ' (', p_cor, ')')) +
        scale_shape_discrete(solid=F, guide=F) + xlab('answer') + ylab(target) +
        scale_colour_discrete(guide=F) + xlim(0, quantile(target_meta$ans_1, 0.95)) + 
        ylim(0, quantile(target_meta$set_count_1, 0.95))
    
    target_meta$ans_1 = target_meta$ans_1 / target_meta$ref_length
    target_meta$set_count_1 = target_meta$set_count_1 / target_meta$set_length
    p_cor = round(cor(target_meta$ans_1, target_meta$set_count_1), 2)
    figure_2 = ggplot(target_meta, aes(ans_1, set_count_1, shape='a', colour="#FF9999"))
    figure_2 = figure_2 + geom_point() + 
        ggtitle(paste0('(normalized) ', organism, ' - ', target_category, ' (', p_cor, ')')) +
        scale_shape_discrete(solid=F, guide=F) + xlab('answer') + ylab(target) +
        scale_colour_discrete(guide=F) + xlim(0, quantile(target_meta$ans_1, 0.95)) + 
        ylim(0, quantile(target_meta$set_count_1, 0.95))
    
    png(paste0('figures/', organism, '_xy_plot_answer_label_', target_category, '_s1.png'), width=1200, height=400)
    multiplot(figure_1, figure_2, cols=2)
    dev.off()
}

# (error) ambiguous
target_meta = subset(meta_data, meta_data$category == 'disaggregate')
target_meta$ans_1 = target_meta$ans_1 / target_meta$ref_length
target_meta$set_count_1 = target_meta$set_count_1 / target_meta$set_length
target_meta$se = sqrt((target_meta$set_count_1 - target_meta$ref_count_1) ^ 2) / (target_meta$ref_count_1 + 1e-10)
cat(quantile(target_meta$se))
figure_1 = ggplot(target_meta, aes(factor(category), se))
figure_1 = figure_1 + geom_boxplot() + 
           ggtitle(paste0(organism, ' - proportional error')) +
           xlab('ambiguous category') + ylab(target) + ylim(0, 5)
png(paste0('figures/', organism, '_error_ambiguous_s1.png'), width=1000, height=400)
print(figure_1)
dev.off()









t = subset(meta_data, meta_data$category != 'disaggregate')

plot(t$ans_1, t$set_count_1)
cor(t$ans_1, t$set_count_1)

plot(meta_data$set_count_1/meta_data$set_length, meta_data$ans_1/meta_data$ref_length)
cor(meta_data$set_count_1, meta_data$ans_1)

path_3 = paste0(organism, '_cdna/', organism, '_trinity_', fold, '_completness.tsv')
completness = read.table(path_3, header=F, stringsAsFactors=F, sep='\t') 
colnames(completness) = c('target_id', 'set_id', 'orientation', 'target_length', 'set_length',
                          'target_alignment', 'set_alignment', 'target_coverage', 'set_coverage')
completness = completness[order(completness$target_id),]
completness$target_count = table(completness$set_id)[completness$set_id]
answer = subset(answer, answer$target_id %in% completness$target_id)
rownames(raw_data_1) = raw_data_1$transcript_id
rownames(raw_data_2) = raw_data_2$transcript_id

path_4 = paste0(organism, '_cdna/blastn_pairwise_', organism, 
               '_cdna_100_ambiguous.tsv')
blast_summary = read.table(path_4, header=T, sep='\t',
                           stringsAsFactors=F)
rownames(blast_summary) =  blast_summary$target_id
completness$target_ambiguous = blast_summary[completness$target_id, 'align_coverage']
    
path_5 = paste0(organism, '_cdna/blastn_pairwise_', organism, 
                '_trinity_10x_ambiguous.tsv')
blast_summary = read.table(path_5, header=T, sep='\t',
                           stringsAsFactors=F)
rownames(blast_summary) =  blast_summary$target_id
completness$set_ambiguous = blast_summary[completness$set_id, 'align_coverage']

#merge_all = cbind(answer[,c(2,3)], raw_data_1$expected_count, raw_data_2$expected_count)
merge_all = cbind(answer[,c(2,3)], raw_data_1[completness$set_id, 'expected_count'], 
                  raw_data_2[completness$set_id, 'expected_count'])
colnames(merge_all) = c('answer_1', 'answer_2', 'rsem_1', 'rsem_2')
rownames(merge_all) = answer$target_id

  
merge_all_t = subset(merge_all, completness$target_count == 1       & 
                                completness$target_ambiguous != 0  &
                                completness$set_ambiguous != 0     &
                                completness$target_coverage > 0.9    &
                                completness$set_coverage > 0.9)

completness$assembly_category = 'others'
condition = completness$target_coverage < 0.2 & completness$set_coverage < 0.2
completness[condition, 'assembly_category'] = 'noise'
condition = completness$target_coverage >= 0.8 & completness$set_coverage >= 0.8
completness[condition, 'assembly_category'] = 'full-length'
condition = completness$target_coverage < 0.2 & completness$set_coverage >= 0.8
completness[condition, 'assembly_category'] = 'fragmented'
condition = completness$target_coverage >= 0.8 & completness$set_coverage < 0.2
completness[condition, 'assembly_category'] = 'over-extended'
condition = completness$target_coverage >= 0.8 & completness$set_coverage >= 0.8
completness[condition, 'assembly_category'] = 'full-length'
condition = completness$target_count > 1
completness[condition, 'assembly_category'] = 'disaggregate'
completness$assembly_category = factor(completness$assembly_category)

print(length(which(completness$assembly_category=='others')))
print(length(which(completness$assembly_category=='noise')))
print(length(which(completness$assembly_category=='full-length')))
print(length(which(completness$assembly_category=='fragmented')))
print(length(which(completness$assembly_category=='over-extended')))
print(length(which(completness$assembly_category=='disaggregate')))

print(length(which(completness$target_ambiguous == 0 & completness$set_ambiguous == 0)))
print(length(which(completness$target_ambiguous != 0 & completness$set_ambiguous == 0)))
print(length(which(completness$target_ambiguous == 0 & completness$set_ambiguous != 0)))
print(length(which(completness$target_ambiguous != 0 & completness$set_ambiguous != 0)))

t = completness[completness$assembly_category=='others', 'target_id']
tt = subset(merge_all, rownames(merge_all) %in% t)
chart.Correlation(tt,  histogram=TRUE, method="pearson", pch=20)

t = completness[which(completness$target_ambiguous != 0 & completness$set_ambiguous == 0), 'target_id']
tt = subset(merge_all, rownames(merge_all) %in% t)
chart.Correlation(tt,  histogram=TRUE, method="pearson", pch=20)

#plot(completness$target_coverage, completness$set_coverage, col=completness$ambiguous_category)

chart.Correlation(merge_all_t, histogram=TRUE, method="pearson", pch=20)

png(filename=paste0('figures/rsem_', organism, '_', fold, '_total.png'), width=800, height=750)  
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
png(filename=paste0('figures/rsem_', organism, '_', fold, '_unique.png'), width=800, height=750)
chart.Correlation(merge_unique, histogram=TRUE, method="pearson", pch=20)
dev.off()

merge_ambiguous = subset(merge_all, rownames(merge_all) %in% ambiguous_set$target_id)
png(filename=paste0('figures/rsem_', organism, '_', fold, '_ambiguous.png'), width=800, height=750)
chart.Correlation(merge_ambiguous, histogram=TRUE, method="pearson", pch=20)
dev.off()

worst = subset(blast_summary, blast_summary$align_coverage > 0.8)
worst_set = subset(merge_all, rownames(merge_all) %in% worst$target_id)
chart.Correlation(worst_set, histogram=TRUE, method="pearson", pch=20)