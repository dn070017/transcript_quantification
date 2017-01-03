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