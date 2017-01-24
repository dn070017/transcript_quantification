library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

###################
# set environment #
###################
{
    setwd('/Users/Hsieh/Documents/Projects/Transcript Quantification/')
    
    organisms = c('yeast', 'mouse')
    folds = c('10x', 'real')
    
    targets = c('reference', 'trinity', 'trinity_answer_label')
    target = targets[2]
    
    prefix = paste0(organism, '_', fold)
    
    out_folder = '20170124'
}

##########################
# set defualt categories #
##########################
{
    assembly_completeness = c('full-length', 'others', 'fragmented', 'over-extended')
    ambiguous_type = c('unique', 'ambiguous')
    ref_set_ambiguous_type = vector(length=4)
    ref_set_nr_ambiguous_type = vector(length=8)
    disaggregate_type = c('single', 'overlap', 'separate')
    
    index = 1
    for(x in seq(1, 2)){
        for(y in seq(1, 2)){
            ref_set_ambiguous_type[index] = paste0(ambiguous_type[x], '_', ambiguous_type[y], sep='')
            index = index + 1
        } 
    }
    
    index = 1
    for(x in seq(1, 2)){
        for(y in seq(1, 2)){
            for(z in seq(1, 2)){
                ref_set_nr_ambiguous_type[index] = paste0(ambiguous_type[x], '_', ambiguous_type[y], '_', ambiguous_type[z], sep='')
                index = index + 1
            }
        } 
    }
}

######################
# create directories #
######################
{
    dir.create(file.path('figures', out_folder), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'all_data'), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'assembly_completeness_type'), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'disaggregate_type'), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'ref_ambiguous_type'), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'set_ambiguous_type'), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'nr_ambiguous_type'), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'ref_set_ambiguous_type'), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'ref_set_nr_ambiguous_type'), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'full_table'), showWarnings=F)
    dir.create(file.path('summary', out_folder), showWarnings=F)
    dir.create(file.path('summary', out_folder, 'all_data'), showWarnings=F)
    dir.create(file.path('summary', out_folder, 'assembly_completeness_type'), showWarnings=F)
    dir.create(file.path('summary', out_folder, 'disaggregate_type'), showWarnings=F)
    dir.create(file.path('summary', out_folder, 'ref_ambiguous_type'), showWarnings=F)
    dir.create(file.path('summary', out_folder, 'set_ambiguous_type'), showWarnings=F)
    dir.create(file.path('summary', out_folder, 'nr_ambiguous_type'), showWarnings=F)
    dir.create(file.path('summary', out_folder, 'ref_set_ambiguous_type'), showWarnings=F)
    dir.create(file.path('summary', out_folder, 'ref_set_nr_ambiguous_type'), showWarnings=F)
    dir.create(file.path('summary', out_folder, 'full_table'), showWarnings=F)
}

######################
# read reference set #
######################
{
    path_1 = paste0(prefix, '/bowtie2_rsem_sample_01_reference/bowtie2_rsem.isoforms.results')
    path_2 = paste0(prefix, '/bowtie2_rsem_sample_02_reference/bowtie2_rsem.isoforms.results')
    raw_data_1 = read.table(path_1, header=T, sep='\t', stringsAsFactors=F)
    raw_data_1 = raw_data_1[order(raw_data_1$transcript_id), c('transcript_id', 'expected_count')]
    raw_data_2 = read.table(path_2, header=T, sep='\t', stringsAsFactors=F)
    raw_data_2 = raw_data_2[order(raw_data_2$transcript_id), c('transcript_id', 'expected_count')]
    rownames(raw_data_1) = raw_data_1$transcript_id
    rownames(raw_data_2) = raw_data_2$transcript_id
    meta_data$ref_count_1 = raw_data_1[meta_data$ref_id, 'expected_count'] + 1
    meta_data$ref_count_2 = raw_data_2[meta_data$ref_id, 'expected_count'] + 1
}

####################
# read exact count #
####################
{
    if(fold != 'real'){
        answer = read.table(paste0(prefix, '/', prefix, '_count_answer.tsv'), header=T, sep='\t', 
                            stringsAsFactors=F)
        answer = answer[order(answer$target_id), ]
        rownames(answer) = answer$target_id
        colnames(answer) = c('target_id', 'answer_1', 'answer_2')
        meta_data$ans_1 = answer[meta_data$ref_id, 'answer_1'] + 1
        meta_data$ans_2 = answer[meta_data$ref_id, 'answer_2'] + 1
    } else {
        meta_data$ans_1 = meta_data$ref_count_1
        meta_data$ans_2 = meta_data$ref_count_2
    }
}

###################    
# read target set #
###################
{
    if(target != 'reference'){
        path_1 = paste0(prefix, '/bowtie2_rsem_sample_01_', target, '/bowtie2_rsem.isoforms.results')
        path_2 = paste0(prefix, '/bowtie2_rsem_sample_02_', target, '/bowtie2_rsem.isoforms.results')
        raw_data_1 = read.table(path_1, header=T, sep='\t', stringsAsFactors=F)
        raw_data_1 = raw_data_1[order(raw_data_1$transcript_id), c('transcript_id', 'expected_count')]
        raw_data_2 = read.table(path_2, header=T, sep='\t', stringsAsFactors=F)
        raw_data_2 = raw_data_2[order(raw_data_2$transcript_id), c('transcript_id', 'expected_count')]
        rownames(raw_data_1) = raw_data_1$transcript_id
        rownames(raw_data_2) = raw_data_2$transcript_id
        meta_data$set_count_1 = raw_data_1[meta_data$set_id, 'expected_count'] + 1
        meta_data$set_count_2 = raw_data_2[meta_data$set_id, 'expected_count'] + 1
    } else {
        meta_data$set_count_1 = meta_data$ref_count_1
        meta_data$set_count_2 = meta_data$ref_count_2
    }
}

###########################
# calculate percent error #
###########################
{
    meta_data$percent_err_1 = round(abs(meta_data$ans_1-meta_data$set_count_1)/meta_data$ans_1, 3)
    meta_data$percent_err_2 = round(abs(meta_data$ans_2-meta_data$set_count_2)/meta_data$ans_2, 3)
}

##################################################
# (figure) quantification scatterplot - all data #
##################################################
{
    p_cor_1 = round(cor(meta_data$ans_1, meta_data$set_count_1), 3)
    s_cor_1 = round(cor(meta_data$ans_1, meta_data$set_count_1, method='spearman'), 3)
    p_cor_2 = round(cor(meta_data$ans_2, meta_data$set_count_2), 3)
    s_cor_2 = round(cor(meta_data$ans_2, meta_data$set_count_2, method='spearman'), 3)
    
    figure_file = paste0('figures/', out_folder, '/all_data/', organism, '_quantification_scatterplot_', target, '_', fold, '.jpg')
    quantification_scatter_plot(meta_data, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1)
    
    stats = data.frame(count=nrow(meta_data), 
                       s1_p_cor=p_cor_1, 
                       s2_p_cor=p_cor_1, 
                       s1_s_cor=s_cor_1,
                       s2_s_cor=s_cor_2,
                       s1_q1_p_err=quantile(meta_data$percent_err_1, 0.25),
                       s1_q2_p_err=quantile(meta_data$percent_err_1, 0.50),
                       s1_q3_p_err=quantile(meta_data$percent_err_1, 0.75),
                       s2_q1_p_err=quantile(meta_data$percent_err_1, 0.25),
                       s2_q2_p_err=quantile(meta_data$percent_err_1, 0.50),
                       s2_q3_p_err=quantile(meta_data$percent_err_1, 0.75))
    summary_file = paste0('summary/', out_folder, '/all_data/', organism, '_', target, '_', fold, '_rsem_summary.tsv')
    write.table(stats, summary_file, sep='\t', quote=F, row.names=F) 
}

####################################################################
# (figure) quantification scatterplot - assembly completeness type #
####################################################################
{
    stats = data.frame()
    
    for(target_category in assembly_completeness) {

        target_meta = subset(meta_data, meta_data$assembly_type == target_category)
        
        if(nrow(target_meta) == 0){
            add = data.frame(assembly_completeness_type=category[[1]][1], 
                             count=0,
                             s1_p_cor=NA, s2_p_cor=NA, s1_s_cor=NA, s2_s_cor=NA,
                             s1_q1_p_err=NA, s1_q2_p_err=NA, s1_q3_p_err=NA,
                             s2_q1_p_err=NA, s2_q2_p_err=NA, s2_q3_p_err=NA)
            if(nrow(stats) == 0){
                stats = add
            } else {
                stats = rbind(stats, add)
            }
            next
        }
        
        p_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1), 3)
        s_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1, method='spearman'), 3)
        p_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2), 3)
        s_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2, method='spearman'), 3)
        
        
        figure_file = paste0('figures/', out_folder, '/assembly_completeness_type/', organism, '_quantification_scatterplot_', target, '_', fold, '_', target_category, '.jpg')
        quantification_scatter_plot(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1)
        
        add = data.frame(assembly_category=target_category, 
                         count=nrow(target_meta), 
                         s1_p_cor=p_cor_1, 
                         s2_p_cor=p_cor_1, 
                         s1_s_cor=s_cor_1,
                         s2_s_cor=s_cor_2,
                         s1_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s1_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s1_q3_p_err=quantile(target_meta$percent_err_1, 0.75),
                         s2_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s2_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s2_q3_p_err=quantile(target_meta$percent_err_1, 0.75))
        if(nrow(stats) == 0){
            stats = add
        } else {
            stats = rbind(stats, add)
        }
    }
    
    summary_file = paste0('summary/', out_folder, '/assembly_completeness_type/', organism, '_', target, '_', fold, '_rsem_summary.tsv')
    write.table(stats, summary_file, sep='\t', quote=F, row.names=F)
}

###########################################################
# (figure) quantification scatterplot - disaggregate type #
###########################################################
{
    stats = data.frame()
    
    for(target_category in disaggregate_type) {

        target_meta = subset(meta_data, meta_data$disaggregate_type == target_category)
        
        if(nrow(target_meta) == 0){
            add = data.frame(disaggregate_type=target_category, 
                             count=0,
                             s1_p_cor=NA, s2_p_cor=NA, s1_s_cor=NA, s2_s_cor=NA,
                             s1_q1_p_err=NA, s1_q2_p_err=NA, s1_q3_p_err=NA,
                             s2_q1_p_err=NA, s2_q2_p_err=NA, s2_q3_p_err=NA)
            if(nrow(stats) == 0){
                stats = add
            } else {
                stats = rbind(stats, add)
            }
            next
        }
        
        p_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1), 3)
        s_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1, method='spearman'), 3)
        p_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2), 3)
        s_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2, method='spearman'), 3)
      
        figure_file = paste0('figures/', out_folder, '/disaggregate_type/', organism, '_quantification_scatterplot_', target, '_', fold, '_', target_category, '.jpg')
        quantification_scatter_plot(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1)
        
        add = data.frame(disaggregate_type=target_category, 
                         count=nrow(target_meta), 
                         s1_p_cor=p_cor_1, 
                         s2_p_cor=p_cor_1, 
                         s1_s_cor=s_cor_1,
                         s2_s_cor=s_cor_2,
                         s1_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s1_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s1_q3_p_err=quantile(target_meta$percent_err_1, 0.75),
                         s2_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s2_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s2_q3_p_err=quantile(target_meta$percent_err_1, 0.75))
        if(nrow(stats) == 0){
            stats = add
        } else {
            stats = rbind(stats, add)
        }
    }
    
    summary_file = paste0('summary/', out_folder, '/disaggregate_type/', organism, '_', target, '_', fold, '_rsem_summary.tsv')
    write.table(stats, summary_file, sep='\t', quote=F, row.names=F)
}

#############################################################
# (figure) quantification scatterplot - reference ambiguous #
#############################################################
{
    stats = data.frame()
    
    for(target_category in ambiguous_type) {
        target_meta = subset(meta_data, meta_data$ref_ambiguous_type == target_category)
        
        if(nrow(target_meta) == 0){
            add = data.frame(ref_ambiguous_type=target_category, 
                             count=0,
                             s1_p_cor=NA, s2_p_cor=NA, s1_s_cor=NA, s2_s_cor=NA,
                             s1_q1_p_err=NA, s1_q2_p_err=NA, s1_q3_p_err=NA,
                             s2_q1_p_err=NA, s2_q2_p_err=NA, s2_q3_p_err=NA)
            if(nrow(stats) == 0){
                stats = add
            } else {
                stats = rbind(stats, add)
            }
            next
        }
        
        p_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1), 3)
        s_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1, method='spearman'), 3)
        p_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2), 3)
        s_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2, method='spearman'), 3)
       
        figure_file = paste0('figures/', out_folder, '/ref_ambiguous_type/', organism, '_quantification_scatterplot_', target, '_', fold, '_', target_category, '.jpg')
        quantification_scatter_plot(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1)
        
        add = data.frame(ref_ambiguous_type=target_category, 
                         count=nrow(target_meta), 
                         s1_p_cor=p_cor_1, 
                         s2_p_cor=p_cor_1, 
                         s1_s_cor=s_cor_1,
                         s2_s_cor=s_cor_2,
                         s1_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s1_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s1_q3_p_err=quantile(target_meta$percent_err_1, 0.75),
                         s2_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s2_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s2_q3_p_err=quantile(target_meta$percent_err_1, 0.75))
        if(nrow(stats) == 0){
            stats = add
        } else {
            stats = rbind(stats, add)
        }
    }
    
    summary_file = paste0('summary/', out_folder, '/ref_ambiguous_type/', organism, '_', target, '_', fold, '_rsem_summary.tsv')
    write.table(stats, summary_file, sep='\t', quote=F, row.names=F)
}

#######################################################
# (figure) quantification scatterplot - set ambiguous #
#######################################################
{
    stats = data.frame()
    
    for(target_category in ambiguous_type) {
        
        target_meta = subset(meta_data, meta_data$set_ambiguous_type == target_category)
        
        if(nrow(target_meta) == 0){
            add = data.frame(set_ambiguous_type=target_category, 
                             count=0,
                             s1_p_cor=NA, s2_p_cor=NA, s1_s_cor=NA, s2_s_cor=NA,
                             s1_q1_p_err=NA, s1_q2_p_err=NA, s1_q3_p_err=NA,
                             s2_q1_p_err=NA, s2_q2_p_err=NA, s2_q3_p_err=NA)
            if(nrow(stats) == 0){
                stats = add
            } else {
                stats = rbind(stats, add)
            }
            next
        }
        
        p_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1), 3)
        s_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1, method='spearman'), 3)
        p_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2), 3)
        s_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2, method='spearman'), 3)
        
        
        figure_file = paste0('figures/', out_folder, '/set_ambiguous_type/', organism, '_quantification_scatterplot_', target, '_', fold, '_', target_category, '.jpg')
        quantification_scatter_plot(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1)
        
        add = data.frame(set_ambiguous_type=target_category, 
                         count=nrow(target_meta), 
                         s1_p_cor=p_cor_1, 
                         s2_p_cor=p_cor_1, 
                         s1_s_cor=s_cor_1,
                         s2_s_cor=s_cor_2,
                         s1_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s1_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s1_q3_p_err=quantile(target_meta$percent_err_1, 0.75),
                         s2_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s2_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s2_q3_p_err=quantile(target_meta$percent_err_1, 0.75))
        if(nrow(stats) == 0){
            stats = add
        } else {
            stats = rbind(stats, add)
        }
    }
    
    summary_file = paste0('summary/', out_folder, '/set_ambiguous_type/', organism, '_', target, '_', fold, '_rsem_summary.tsv')
    write.table(stats, summary_file, sep='\t', quote=F, row.names=F)
}

#######################################################
# (figure) quantification scatterplot - nr ambiguous #
#######################################################
{
    stats = data.frame()
    
    for(target_category in ambiguous_type) {

        target_meta = subset(meta_data, meta_data$nr_ambiguous_type == target_category)
        
        if(nrow(target_meta) == 0){
            add = data.frame(nr_ambiguous_type=target_category, 
                             count=0,
                             s1_p_cor=NA, s2_p_cor=NA, s1_s_cor=NA, s2_s_cor=NA,
                             s1_q1_p_err=NA, s1_q2_p_err=NA, s1_q3_p_err=NA,
                             s2_q1_p_err=NA, s2_q2_p_err=NA, s2_q3_p_err=NA)
            if(nrow(stats) == 0){
                stats = add
            } else {
                stats = rbind(stats, add)
            }
            next
        }
        
        p_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1), 3)
        s_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1, method='spearman'), 3)
        p_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2), 3)
        s_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2, method='spearman'), 3)
        
        
        figure_file = paste0('figures/', out_folder, '/nr_ambiguous_type/', organism, '_quantification_scatterplot_', target, '_', fold, '_', target_category, '.jpg')
        quantification_scatter_plot(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1)
        
        add = data.frame(nr_ambiguous_type=target_category, 
                         count=nrow(target_meta), 
                         s1_p_cor=p_cor_1, 
                         s2_p_cor=p_cor_1, 
                         s1_s_cor=s_cor_1,
                         s2_s_cor=s_cor_2,
                         s1_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s1_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s1_q3_p_err=quantile(target_meta$percent_err_1, 0.75),
                         s2_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s2_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s2_q3_p_err=quantile(target_meta$percent_err_1, 0.75))
        if(nrow(stats) == 0){
            stats = add
        } else {
            stats = rbind(stats, add)
        }
    }
    
    summary_file = paste0('summary/', out_folder, '/nr_ambiguous_type/', organism, '_', target, '_', fold, '_rsem_summary.tsv')
    write.table(stats, summary_file, sep='\t', quote=F, row.names=F)
}

##################################################################
# (figure) quantification scatterplot - reference, set ambiguous #
##################################################################
{
    meta_data$category = with(answer_label, paste(ref_ambiguous_type, set_ambiguous_type, 
                                                  sep='_'))
    stats = data.frame()
    
    for(target_category in ref_set_ambiguous_type) {

        target_meta = subset(meta_data, meta_data$category == target_category)
        
        if(nrow(target_meta) == 0){
            category = strsplit(target_category, '_')
            add = data.frame(ref_ambiguous_type=category[[1]][1],
                             set_ambiguous_type=category[[1]][2],
                             count=0,
                             s1_p_cor=NA, s2_p_cor=NA, s1_s_cor=NA, s2_s_cor=NA,
                             s1_q1_p_err=NA, s1_q2_p_err=NA, s1_q3_p_err=NA,
                             s2_q1_p_err=NA, s2_q2_p_err=NA, s2_q3_p_err=NA)
            if(nrow(stats) == 0){
                stats = add
            } else {
                stats = rbind(stats, add)
            }
            next
        }
        
        p_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1), 3)
        s_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1, method='spearman'), 3)
        p_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2), 3)
        s_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2, method='spearman'), 3)
        
        
        figure_file = paste0('figures/', out_folder, '/ref_set_ambiguous_type/', organism, '_quantification_scatterplot_', target, '_', fold, '_', target_category, '.jpg')
        quantification_scatter_plot(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1)
        
        category = strsplit(target_category, '_')
        add = data.frame(ref_ambiguous_type=category[[1]][1],
                         set_ambiguous_type=category[[1]][2],
                         count=nrow(target_meta), 
                         s1_p_cor=p_cor_1, 
                         s2_p_cor=p_cor_1, 
                         s1_s_cor=s_cor_1,
                         s2_s_cor=s_cor_2,
                         s1_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s1_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s1_q3_p_err=quantile(target_meta$percent_err_1, 0.75),
                         s2_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s2_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s2_q3_p_err=quantile(target_meta$percent_err_1, 0.75))
        if(nrow(stats) == 0){
            stats = add
        } else {
            stats = rbind(stats, add)
        }
    }
    
    summary_file = paste0('summary/', out_folder, '/ref_set_ambiguous_type/', organism, '_', target, '_', fold, '_rsem_summary.tsv')
    write.table(stats, summary_file, sep='\t', quote=F, row.names=F)
}

######################################################################
# (figure) quantification scatterplot - reference, set, nr ambiguous #
######################################################################
{
    meta_data$category = with(answer_label, paste(ref_ambiguous_type, set_ambiguous_type, 
                                                  nr_ambiguous_type, sep='_'))
    stats = data.frame()
    
    for(target_category in ref_set_nr_ambiguous_type) {

        target_meta = subset(meta_data, meta_data$category == target_category)
        
        if(nrow(target_meta) == 0){
            category = strsplit(target_category, '_')
            add = data.frame(ref_ambiguous_type=category[[1]][1],
                             set_ambiguous_type=category[[1]][2],
                             nr_ambiguous_type=category[[1]][3],
                             count=0,
                             s1_p_cor=NA, s2_p_cor=NA, s1_s_cor=NA, s2_s_cor=NA,
                             s1_q1_p_err=NA, s1_q2_p_err=NA, s1_q3_p_err=NA,
                             s2_q1_p_err=NA, s2_q2_p_err=NA, s2_q3_p_err=NA)
            if(nrow(stats) == 0){
                stats = add
            } else {
                stats = rbind(stats, add)
            }
            next
        }
        
        p_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1), 3)
        s_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1, method='spearman'), 3)
        p_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2), 3)
        s_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2, method='spearman'), 3)
        
        
        figure_file = paste0('figures/', out_folder, '/ref_set_nr_ambiguous_type/', organism, '_quantification_scatterplot_', target, '_', fold, '_', target_category, '.jpg')
        quantification_scatter_plot(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1)
        
        category = strsplit(target_category, '_')
        add = data.frame(ref_ambiguous_type=category[[1]][1],
                         set_ambiguous_type=category[[1]][2],
                         nr_ambiguous_type=category[[1]][3],
                         count=nrow(target_meta), 
                         s1_p_cor=p_cor_1, 
                         s2_p_cor=p_cor_1, 
                         s1_s_cor=s_cor_1,
                         s2_s_cor=s_cor_2,
                         s1_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s1_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s1_q3_p_err=quantile(target_meta$percent_err_1, 0.75),
                         s2_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s2_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s2_q3_p_err=quantile(target_meta$percent_err_1, 0.75))
        if(nrow(stats) == 0){
            stats = add
        } else {
            stats = rbind(stats, add)
        }
    }
    
    summary_file = paste0('summary/', out_folder, '/ref_set_nr_ambiguous_type/', organism, '_', target, '_', fold, '_rsem_summary.tsv')
    write.table(stats, summary_file, sep='\t', quote=F, row.names=F)
}

####################################################
# (figure) quantification scatterplot - full table #
####################################################
{
    meta_data$category = with(answer_label, paste(assembly_type, disaggregate_type, 
                                                  ref_ambiguous_type, set_ambiguous_type, 
                                                  nr_ambiguous_type, sep='_'))
    stats = data.frame()
    
    for(target_category in unique(meta_data$category)) {
        
        target_meta = subset(meta_data, meta_data$category == target_category)
        
        if(nrow(target_meta) == 0){
            category = strsplit(target_category, '_')
            add = data.frame(assembly_category=category[[1]][1], 
                             disaggregate_type=category[[1]][2], 
                             ref_ambiguous_type=category[[1]][3],
                             set_ambiguous_type=category[[1]][4],
                             nr_ambiguous_type=category[[1]][5],
                             count=0,
                             s1_p_cor=NA, s2_p_cor=NA, s1_s_cor=NA, s2_s_cor=NA,
                             s1_q1_p_err=NA, s1_q2_p_err=NA, s1_q3_p_err=NA,
                             s2_q1_p_err=NA, s2_q2_p_err=NA, s2_q3_p_err=NA)
            if(nrow(stats) == 0){
                stats = add
            } else {
                stats = rbind(stats, add)
            }
            next
        }
        
        p_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1), 3)
        s_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1, method='spearman'), 3)
        p_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2), 3)
        s_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2, method='spearman'), 3)
        
        
        figure_file = paste0('figures/', out_folder, '/full_table/', organism, '_quantification_scatterplot_', target, '_', fold, '_', target_category, '.jpg')
        quantification_scatter_plot(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1)
        
        category = strsplit(target_category, '_')
        add = data.frame(assembly_category=category[[1]][1], 
                         disaggregate_type=category[[1]][2], 
                         ref_ambiguous_type=category[[1]][3],
                         set_ambiguous_type=category[[1]][4],
                         nr_ambiguous_type=category[[1]][5],
                         count=nrow(target_meta), 
                         s1_p_cor=p_cor_1, 
                         s2_p_cor=p_cor_1, 
                         s1_s_cor=s_cor_1,
                         s2_s_cor=s_cor_2,
                         s1_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s1_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s1_q3_p_err=quantile(target_meta$percent_err_1, 0.75),
                         s2_q1_p_err=quantile(target_meta$percent_err_1, 0.25),
                         s2_q2_p_err=quantile(target_meta$percent_err_1, 0.50),
                         s2_q3_p_err=quantile(target_meta$percent_err_1, 0.75))
        if(nrow(stats) == 0){
            stats = add
        } else {
            stats = rbind(stats, add)
        }
    }
    
    summary_file = paste0('summary/', out_folder, '/full_table/', organism, '_', target, '_', fold, '_rsem_summary.tsv')
    write.table(stats, summary_file, sep='\t', quote=F, row.names=F)
}

#################################
# percent error cumulative plot #
#################################
{
    component_set = c('all_data', 'assembly_completeness_type', 'disaggregate_type', 
                      'ref_ambiguous_type', 'set_ambiguous_type', 'nr_ambiguous_type',
                      'ref_set_ambiguous_type', 'ref_set_nr_ambiguous_type')
    for(component in component_set){
        figure_file = paste0('figures/', out_folder, '/', component, '/', organism, '_percent_error_cumulative_plot_', target, '_', fold, '.jpg')
        if(component == 'all_data'){
            figure_1 = ggplot(meta_data, aes(percent_err_1, shape='a'))
            figure_2 = ggplot(meta_data, aes(percent_err_2, shape='a'))
        } else if (component == 'assembly_completeness_type'){
            figure_1 = ggplot(meta_data, aes(percent_err_1, shape='a', colour=assembly_type))
            figure_2 = ggplot(meta_data, aes(percent_err_2, shape='a', colour=assembly_type))
        } else if (component == 'disaggregate_type'){
            figure_1 = ggplot(meta_data, aes(percent_err_1, shape='a', colour=disaggregate_type))
            figure_2 = ggplot(meta_data, aes(percent_err_2, shape='a', colour=disaggregate_type))
        } else if (component == 'ref_ambiguous_type') {
            figure_1 = ggplot(meta_data, aes(percent_err_1, shape='a', colour=ref_ambiguous_type))
            figure_2 = ggplot(meta_data, aes(percent_err_2, shape='a', colour=ref_ambiguous_type))
        } else if (component == 'set_ambiguous_type') {
            figure_1 = ggplot(meta_data, aes(percent_err_1, shape='a', colour=set_ambiguous_type))
            figure_2 = ggplot(meta_data, aes(percent_err_2, shape='a', colour=set_ambiguous_type))
        } else if (component == 'nr_ambiguous_type') {
            figure_1 = ggplot(meta_data, aes(percent_err_1, shape='a', colour=nr_ambiguous_type))
            figure_2 = ggplot(meta_data, aes(percent_err_2, shape='a', colour=nr_ambiguous_type))
        } else if (component == 'ref_set_ambiguous_type') {
            meta_data$category = with(answer_label, paste(ref_ambiguous_type, set_ambiguous_type, sep='_'))
            figure_1 = ggplot(meta_data, aes(percent_err_1, shape='a', colour=category))
            figure_2 = ggplot(meta_data, aes(percent_err_2, shape='a', colour=category))
        } else if (component == 'ref_set_nr_ambiguous_type') {
            meta_data$category = with(answer_label, paste(ref_ambiguous_type, set_ambiguous_type, nr_ambiguous_type, sep='_'))
            figure_1 = ggplot(meta_data, aes(percent_err_1, shape='a', colour=category))
            figure_2 = ggplot(meta_data, aes(percent_err_2, shape='a', colour=category))
        }
        figure_1 = figure_1 + stat_ecdf() + 
                   scale_shape_discrete(solid=F, guide=F) + 
                   xlab('percent error threshold') + 
                   ylab('% of data smaller than percent error threshold') +
                   scale_x_continuous(expand=c(0.005, 0), limits=c(0, 1)) +
                   scale_y_continuous(expand=c(0.005, 0), limits=c(0, 1))
        
        figure_2 = figure_2 + stat_ecdf() + 
                   scale_shape_discrete(solid=F, guide=F) + 
                   xlab('percent error threshold') + 
                   ylab('% of data smaller than percent error threshold') +
                   scale_x_continuous(expand=c(0.005, 0), limits=c(0, 1)) +
                   scale_y_continuous(expand=c(0.005, 0), limits=c(0, 1))
        if(component == 'all_data'){
            jpeg(figure_file, res=350, width=3000, height=1500)
            multiplot(figure_1, figure_2, cols=2)
            dev.off()
        } else {
            jpeg(figure_file, res=350, width=3000, height=1500)
            grid_arrange_shared_legend(figure_1, figure_2, ncol=2, nrow=1)
            dev.off()
        }
    }
}

###################################################
# (archive) scatterplot (percet error - category) #
###################################################
{
    # p_cor_3 = round(cor(meta_data$set_count, meta_data$p_err_1), 3)
    # s_cor_3 = round(cor(meta_data$set_count, meta_data$p_err_1, method='spearman'), 3)
    # p_cor_4 = round(cor(meta_data$set_count, meta_data$p_err_2), 3)
    # s_cor_4 = round(cor(meta_data$set_count, meta_data$p_err_2, method='spearman'), 3)
    # 
    # figure_3 = ggplot(meta_data, aes(set_count, p_err_1, shape='a', colour="#FF9999"))
    # figure_3 = figure_3 + geom_point() + 
    #     annotate("text", x=0, y=Inf, hjust = 0, vjust=1.25, size=3.5,
    #              label=paste0("Pearson's R: ", p_cor_3, '\n', "Spearman's R: ", s_cor_3, sep='')) +
    #     scale_shape_discrete(solid=F, guide=F) + 
    #     #xlab('ambiguously mapping coverage') + 
    #     #xlab('assembly coverage') +
    #     xlab('disaggregate transcripts number') + 
    #     ylab('percent error') +
    #     scale_colour_discrete(guide=F) + 
    #     #xlim(0, 1) + 
    #     ylim(0, 1)
    # 
    # figure_4 = ggplot(meta_data, aes(set_count, p_err_2, shape='a', colour="#FF9999"))
    # figure_4 = figure_4 + geom_point() + 
    #     annotate("text", x=0, y=Inf, hjust = 0, vjust=1.25, size=3.5,
    #              label=paste0("Pearson's R: ", p_cor_4, '\n', "Spearman's R: ", s_cor_4, sep='')) +
    #     scale_shape_discrete(solid=F, guide=F) + 
    #     #xlab('ambiguously mapping coverage') +
    #     #xlab('assembly coverage') +
    #     xlab('disaggregate transcripts number') +         
    #     ylab('percent error') +
    #     scale_colour_discrete(guide=F) + 
    #     #xlim(0, 1) +
    #     ylim(0, 1)
    # 
    # png(paste0('figures/', out_folder, '/', organism, '_perr_disaggregate_', target, '_', fold, '.png'), res=350, width=3000, height=1500)
    # multiplot(figure_3, figure_4, cols=2)
    # dev.off()
}
