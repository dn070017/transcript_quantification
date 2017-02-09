library('ggplot2')
library('randomForest')
library('reprtree')

###################
# set environment #
###################
{
    setwd('/Users/Hsieh/Documents/Projects/Transcript Quantification/')
    
    organisms <- c('yeast', 'mouse')
    organism <- organisms[2]
    
    targets <- c('trinity', 'split')
    target <- targets[1]
    
    folds <- c('10x', 'real')
    fold <- folds[1]
    
    out_folder = '20170124'
    dir.create(file.path('figures', out_folder), showWarnings=F)
    dir.create(file.path('figures', out_folder, 'manuscript'), showWarnings=F)
}

############################
# part 1: reference branch #
############################
{
    target_meta = meta_data[, c('percent_err_1', 'percent_err_2', 'ref_ambiguous_type')]
    colnames(target_meta) = c('percent_err_1', 'percent_err_2', 'category')
    figure_1 = ggplot(target_meta, aes(percent_err_1, shape='a', colour=category))
    figure_1 = figure_1 + stat_ecdf() + 
               scale_shape_discrete(solid=F, guide=F) + 
               ggtitle('A') +
               xlab('percentage error threshold') + 
               ylab('proportion of transcripts with an\nerror smaller than the threshold') +
               scale_x_continuous(expand=c(0.010, 0), limits=c(0, 1)) +
               scale_y_continuous(expand=c(0.005, 0), limits=c(0, 1))
    
    figure_2 = ggplot(target_meta, aes(percent_err_2, shape='a', colour=category))
    figure_2 = figure_2 + stat_ecdf() + 
               scale_shape_discrete(solid=F, guide=F) + 
               ggtitle('B') +
               xlab('percentage error threshold') + 
               ylab('proportion of transcripts with an\nerror smaller than the threshold') +
               scale_x_continuous(expand=c(0.010, 0), limits=c(0, 1)) +
               scale_y_continuous(expand=c(0.005, 0), limits=c(0, 1))
    figure_file = paste0('figures/', out_folder, '/manuscript/01_', organism, '_pec_ra_', target, '_', fold, '.jpg')
    jpeg(figure_file, res=350, width=3000, height=1200)
    grid_arrange_shared_legend(figure_1, figure_2, ncol=2, nrow=1)
    dev.off()
    
    for(target_category in ambiguous_type) {
        target_meta = subset(meta_data, meta_data$ref_ambiguous_type == target_category)
        
        p_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1), 3)
        s_cor_1 = round(cor(target_meta$ans_1, target_meta$set_count_1, method='spearman'), 3)
        p_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2), 3)
        s_cor_2 = round(cor(target_meta$ans_2, target_meta$set_count_2, method='spearman'), 3)
        
        figure_file = paste0('figures/', out_folder, '/manuscript/01_', organism, '_qs_ra_', target, '_', fold, '_', target_category, '.jpg')
        if(target_category == 'unique'){
            quantification_scatter_plot_manuscript(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1, 'C', 'D')
        } else {
            quantification_scatter_plot_manuscript(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_1, s_cor_1, 'E', 'F')
        } 
    } 
}

#################################
# part 3: assembly completeness #
#################################
# A. assembly completeness #
{
    # re-run read_meta_data.R first #
    target_meta = answer_label
    assembly_type_count = table(target_meta$assembly_type)
    
    all = nrow(target_meta)
    fr = paste0('fragmented (', round(assembly_type_count['fragmented']/all, 3), ')')
    fu = paste0('full-length (', round(assembly_type_count['full-length']/all, 3), ')')
    no = paste0('noise (', round(assembly_type_count['noise']/all, 3), ')')
    oe = paste0('over-extended (', round(assembly_type_count['over-extended']/all, 3), ')')
    ot = paste0('others (', round(assembly_type_count['others']/all, 3), ')')
    levels(target_meta$assembly_type) = c(fr, fu, no, ot, oe)
    target_meta$category = target_meta$assembly_type
    figure = ggplot(target_meta, aes(ref_coverage, set_coverage, shape='a', colour=category))
    figure = figure + geom_point() + 
             scale_shape_discrete(solid=F, guide=F) + 
             expand_limits(x=0, y=0) +
             #guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
             #ggtitle('A') +
             #theme(legend.position="bottom",legend.direction="horizontal") +
             scale_x_continuous(expand=c(0.025, 0), limits=c(0, 1)) +
             scale_y_continuous(expand=c(0.025, 0), limits=c(0, 1)) +
             xlab('transcript coverage') +
             ylab('contig coverage')
     
    figure_file = paste0('figures/', out_folder, '/manuscript/02_ac_all_', organism, '_', target, '_', fold, '.jpg')
    jpeg(figure_file, res=350, width=2000, height=1200)
    multiplot(figure, cols=1)
    dev.off()
    
    # target_meta = subset(answer_label, answer_label$ref_ambiguous_site == 0)
    # assembly_type_count = table(target_meta$assembly_type)
    # all = nrow(target_meta)
    # fr = paste0('fragmented (', round(assembly_type_count['fragmented']/all, 3), ')')
    # fu = paste0('full-length (', round(assembly_type_count['full-length']/all, 3), ')')
    # no = paste0('noise (', round(assembly_type_count['noise']/all, 3), ')')
    # oe = paste0('over-extended (', round(assembly_type_count['over-extended']/all, 3), ')')
    # ot = paste0('others (', round(assembly_type_count['others']/all, 3), ')')
    # levels(target_meta$assembly_type) = c(fr, fu, no, ot, oe)
    # figure = ggplot(target_meta, aes(ref_coverage, set_coverage, shape='a', colour=assembly_type))
    # figure = figure + geom_point() + 
    #          scale_shape_discrete(solid=F, guide=F) + 
    #          expand_limits(x=0, y=0) +
    #          guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    #          ggtitle('B') +
    #          theme(legend.position="bottom",legend.direction="horizontal") +
    #          scale_x_continuous(expand=c(0.025, 0), limits=c(0, 1)) +
    #          scale_y_continuous(expand=c(0.025, 0), limits=c(0, 1)) +
    #          xlab('reference coverage') +
    #          ylab('assembly coverage')
    # 
    # figure_file = paste0('figures/', out_folder, '/manuscript/02_ac_ref_uniq_', organism, '_', target, '_', fold, '.jpg')
    # jpeg(figure_file, res=350, width=2000, height=1500)
    # multiplot(figure, cols=1)
    # dev.off()
    # 
    # target_meta = subset(answer_label, answer_label$ref_ambiguous_site != 0)
    # assembly_type_count = table(target_meta$assembly_type)
    # all = nrow(target_meta)
    # fr = paste0('fragmented (', round(assembly_type_count['fragmented']/all, 3), ')')
    # fu = paste0('full-length (', round(assembly_type_count['full-length']/all, 3), ')')
    # no = paste0('noise (', round(assembly_type_count['noise']/all, 3), ')')
    # oe = paste0('over-extended (', round(assembly_type_count['over-extended']/all, 3), ')')
    # ot = paste0('others (', round(assembly_type_count['others']/all, 3), ')')
    # levels(target_meta$assembly_type) = c(fr, fu, no, ot, oe)
    # figure = ggplot(target_meta, aes(ref_coverage, set_coverage, shape='a', colour=assembly_type))
    # figure = figure + geom_point() + 
    #          scale_shape_discrete(solid=F, guide=F) + 
    #          expand_limits(x=0, y=0) +
    #          guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    #          ggtitle('C') +
    #          theme(legend.position="bottom",legend.direction="horizontal") +
    #          scale_x_continuous(expand=c(0.025, 0), limits=c(0, 1)) +
    #          scale_y_continuous(expand=c(0.025, 0), limits=c(0, 1)) +
    #          xlab('reference coverage') +
    #          ylab('assembly coverage')
    # 
    # figure_file = paste0('figures/', out_folder, '/manuscript/02_ac_ref_amb_', organism, '_', target, '_', fold, '.jpg')
    # jpeg(figure_file, res=350, width=2000, height=1500)
    # multiplot(figure, cols=1)
    # dev.off()
}
# B. assembly completness error cumulative #
{
    library(RColorBrewer)
    meta_data$percent_err = (meta_data$percent_err_1 + meta_data$percent_err_2) / 2 
    meta_data$category_am = with(meta_data, paste(ref_ambiguous_type, set_ambiguous_type, sep='_'))
    meta_data$category = meta_data$assembly_type
    index = 1
    figures = vector("list", 4) 
    my_colors = c('#FF9999', '#7AAC45', '#26ED29', '#04F7FE', '#FA04FE')
    names(my_colors) = c('fragmented', 'full-length', 'noise', 'others', 'over-extended')
    col_scale = scale_color_manual(name='category', values=my_colors)
    category_types = c('ambiguous_unique', 'unique_unique', 'ambiguous_ambiguous', 'unique_ambiguous')
    #category_types = c('unique_unique', 'unique_ambiguous', 'ambiguous_unique', 'ambiguous_ambiguous')
    for(target_category in category_types){
        id = 'A'
        if(index == 1) { id = 'A' }
        if(index == 2) { id = 'B' }
        if(index == 3) { id = 'C' }
        if(index == 4) { id = 'D' }
        target_meta = meta_data
        target_meta = subset(target_meta, target_meta$category_am == target_category)

        figure = ggplot(target_meta, aes(percent_err, shape='a', colour=category))
        figure = figure + stat_ecdf() + 
                 scale_shape_discrete(solid=F, guide=F) + 
                 col_scale +
                 ggtitle(id) +
                 xlab('percentage error threshold') + 
                 ylab('proportion of transcripts with an\nerrorsmaller than the threshold') +
                 scale_x_continuous(expand=c(0.010, 0), limits=c(0, 1)) +
                 scale_y_continuous(expand=c(0.005, 0), limits=c(0, 1))
        figures[[index]] = figure
        index = index + 1
    }
    figure_file = paste0('figures/', out_folder, '/manuscript/02_ac_all_', organism, '_pec_ra_', target, '_', fold, '.jpg')
    jpeg(figure_file, res=350, width=2800, height=2800)
    grid_arrange_shared_legend(figures[[1]], figures[[2]], figures[[3]], figures[[4]], ncol=2, nrow=2)
    dev.off()
}

###############################
# part 4: disaggreagate types #
###############################
{
    meta_data$category = meta_data$disaggregate_type
    meta_data$percent_err = (meta_data$percent_err_1 + meta_data$percent_err_2) / 2
    figure = ggplot(meta_data, aes(percent_err, shape='a', colour=category))
    figure = figure + stat_ecdf() + 
             scale_shape_discrete(solid=F, guide=F) + 
             #ggtitle('A') +
             xlab('percentage error threshold') + 
             ylab('proportion of transcripts with an\nerror smaller than the threshold') +
             scale_x_continuous(expand=c(0.010, 0), limits=c(0, 1)) +
             scale_y_continuous(expand=c(0.005, 0), limits=c(0, 1))
    figure_file = paste0('figures/', out_folder, '/manuscript/03_dt_all_', organism, '_pec_ra_', target, '_', fold, '.jpg')
    jpeg(figure_file, res=350, width=1500, height=1200)
    grid_arrange_shared_legend(figure, ncol=1)
    dev.off()  
}

################################################################
# chi-square test for ambiguous site and assembly completeness #
################################################################
{
    library(MASS)
    target_meta = meta_data
    target_meta$category = with(target_meta, paste(ref_ambiguous_type, set_ambiguous_type, sep='_'))
    table_assembly_ambiguos = table(target_meta$category, target_meta$assembly_type)
    chisq.test(table_assembly_ambiguos) 
}

######################################################
# chi-square test for assembly type and disaggregate #
######################################################
{
    library(MASS)
    target_meta = meta_data
    target_assembly_disaggregate = table(target_meta$assembly_type, target_meta$disaggregate_type)
    target_assembly_disaggregate = target_assembly_disaggregate[c(1,2,4,5), ]
    chisq.test(target_assembly_disaggregate)
}

################################
# Principle Component Analysis #
################################
{
    meta_data$percent_err = (meta_data$percent_err_1 + meta_data$percent_err_2) / 2
    features = c(4, 5, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20)   
    target_meta = meta_data[, features]
    colnames(target_meta) = c('ref_len', 'contig_len', 'ref_cov', 'contig_cov',
                              '#_same_contig', 'ref_amb_site', 'ref_amb_cov', 'contig_amb_site',
                              'contig_amb_cov', 'single', 'separate', 'overlap')
    target_pca = prcomp(target_meta, center = TRUE, scale. = TRUE)
    df = data.frame(x=target_pca$rotation[,1], y=target_pca$rotation[, 2], group=rownames(target_pca$rotation))
    
    figure = ggplot(df, aes(x, y, shape='a', colour=group, label=group))
    figure = figure + geom_point() + geom_text(nudge_y=0.03) +
             scale_shape_discrete(solid=F, guide=F) + 
             scale_colour_discrete(guide=F) +
             #ggtitle('A') +
             xlim(-0.6, 0.6) +
             xlab('PC1') + 
             ylab('PC2')
             #scale_x_continuous(expand=c(0.010, 0), limits=c(0, 1)) +
             #scale_y_continuous(expand=c(0.005, 0), limits=c(0, 1))
    figure_file = paste0('figures/', out_folder, '/manuscript/04_pca_', organism, '_', target, '_', fold, '.jpg')
    jpeg(figure_file, res=350, width=2000, height=2000)
    print(figure, ncol=1)
    dev.off()  
}
# archive: random forest #
##########################
{
    # features = c(4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20)
    # 
    # analyze = meta_data[,features]
    # 
    # #for(x in seq(1, 14)){
    # #    analyze[,x] = (analyze[, x] - min(analyze[, x])) / (max(analyze[,x]) - min(analyze[,x]))
    # #}
    # 
    # error = (meta_data$percent_err_1 + meta_data$percent_err_2) / 2
    # # analyze$error = error
    # 
    # analyze$error_type = 's'
    # condition = error > 0.10 & error <= 0.25
    # analyze[condition, 'error_type'] = 'a'
    # condition = error > 0.25 & error <= 0.5
    # analyze[condition, 'error_type'] = 'b'
    # condition = error > 0.5 & error <= 0.75
    # analyze[condition, 'error_type'] = 'c'
    # condition = error > 0.75
    # analyze[condition, 'error_type'] = 'f'
    # analyze$error_type = factor(analyze$error_type)
    # 
    # analyze = analyze[1:1000,]
    # rf = randomForest(error_type ~ ., data=analyze, ntree=5, proximity=T,
    #                   keep.forest=TRUE, importance=TRUE, mtry=2, do.trace=100)
    # png('test.png', width=15000, height=5000, res=300)
    # reprtree:::plot.getTree(rf)
    # dev.off()
}
