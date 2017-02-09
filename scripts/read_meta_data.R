library(ggplot2)

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
    dir.create(file.path('figures', out_folder, 'assembly_completeness_type'), showWarnings=F)
}

##############################
# read assembly completeness #
##############################
{
    if (target == 'trinity') {
        answer_label_path <- paste0(organism, '_data/answer_label_', organism, '_trinity_', fold, '.tsv')
    } else if (target == 'split') {
        answer_label_path <- paste0(organism, '_data/answer_label_', organism, '_trinity_', fold, '_answer_label_split.tsv')
    }
    answer_label <- read_answer_label(answer_label_path)
}

#######################
# read ambiguous site #
#######################
{
    # 1. reference #
    ambiguous_ref_path = paste0(organism, '_data/ambiguous_site_', organism, '_cdna.tsv')
    answer_label <- read_ambiguous_sites(ambiguous_ref_path, answer_label, 'ref')
    # 2. trinity and nr #
    if (target == 'trinity') {
        ambiguous_set_path = paste0(organism, '_data/ambiguous_site_', organism, '_trinity_', fold, '.tsv')
        ambiguous_nr_path = paste0(organism, '_data/ambiguous_site_', organism, '_trinity_', fold, '_answer_label.tsv')
        answer_label <- read_ambiguous_sites(ambiguous_set_path, answer_label, 'set')
        answer_label <- read_ambiguous_sites(ambiguous_nr_path, answer_label, 'nr')
    } else if (target == 'split') {
        ambiguous_set_path = paste0(organism, '_data/ambiguous_site_', organism, '_trinity_', fold, '_answer_label_split.tsv')
        ambiguous_nr_path = paste0(organism, '_data/ambiguous_site_', organism, '_trinity_', fold, '_answer_label_split_nr.tsv')
        answer_label <- read_ambiguous_sites(ambiguous_set_path, answer_label, 'set')
        answer_label <- read_ambiguous_sites(ambiguous_nr_path, answer_label, 'nr')
    }
}

##########################
# classify assembly type #
##########################
{
    answer_label$assembly_type = 'others'
    condition = answer_label$ref_coverage < 0.5 & answer_label$set_coverage < 0.25
    answer_label[condition, 'assembly_type'] = 'noise'
    condition = answer_label$ref_coverage < 0.25 & answer_label$set_coverage < 0.5
    answer_label[condition, 'assembly_type'] = 'noise'
    condition = answer_label$ref_coverage >= 0.75 & answer_label$set_coverage >= 0.75
    answer_label[condition, 'assembly_type'] = 'full-length'
    condition = answer_label$ref_coverage < 0.25 & answer_label$set_coverage >= 0.5
    answer_label[condition, 'assembly_type'] = 'fragmented'
    condition = answer_label$ref_coverage >= 0.5 & answer_label$set_coverage < 0.25
    answer_label[condition, 'assembly_type'] = 'over-extended'
    answer_label$assembly_type = factor(answer_label$assembly_type)
}

########################################
# (figures) assembly completeness type #
########################################
{
    figure = ggplot(answer_label, aes(ref_coverage, set_coverage, shape='a', colour=assembly_type))
    figure = figure + geom_point() + 
             scale_shape_discrete(solid=F, guide=F) + 
             expand_limits(x=0, y=0) +
             scale_x_continuous(expand=c(0.025, 0), limits=c(0, 1)) +
             scale_y_continuous(expand=c(0.025, 0), limits=c(0, 1)) +
             xlab('reference coverage') +
             ylab('assembly coverage')
             
    figure_file = paste0('figures/', out_folder, '/assembly_completeness_type/assembly_completeness_type_', organism, '_', target, '_', fold, '.jpg')
    jpeg(figure_file, res=350, width=2000, height=1500)
    multiplot(figure, cols=1)
    dev.off()
}

#####################
# remove noise data #
#####################
answer_label = subset(answer_label, answer_label$assembly_type != 'noise')

##########################
# read disaggregate type #
##########################
{
    if (target == 'trinity') {
        disaggregate_type_path = paste0(organism, '_data/disaggregate_type_', organism, '_trinity_', fold, '.tsv')
        answer_label = read_disaggregate_type(disaggregate_type_path, answer_label)
    } else if (target == 'split') {
        disaggregate_type_path = paste0(organism, '_data/disaggregate_type_', organism, '_trinity_', fold, '_answer_label_split.tsv')
        answer_label = read_disaggregate_type(disaggregate_type_path, answer_label)
    }
}

##############################
# classify disaggregate type #
##############################
{
    answer_label$disaggregate_type = 'overlap'
    condition = answer_label$disaggregate_type_unique > 0
    answer_label[condition, 'disaggregate_type'] = 'single'
    condition = answer_label$disaggregate_type_separate > 0
    answer_label[condition, 'disaggregate_type'] = 'separate'
}

###########################
# classify ambiguous type #
###########################
{
    answer_label$ref_ambiguous_type = 'unique'
    condition = answer_label$ref_ambiguous_site > 0
    answer_label[condition, 'ref_ambiguous_type'] = 'ambiguous'
    answer_label$set_ambiguous_type = 'unique'
    condition = answer_label$set_ambiguous_site > 0
    answer_label[condition, 'set_ambiguous_type'] = 'ambiguous'
    answer_label$nr_ambiguous_type = 'unique'
    condition = answer_label$nr_ambiguous_site > 0
    answer_label[condition, 'nr_ambiguous_type'] = 'ambiguous'
}

meta_data = answer_label
