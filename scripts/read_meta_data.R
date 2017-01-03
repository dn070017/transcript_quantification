library(ggplot2)

setwd('/Users/Hsieh/Documents/Projects/Transcript Quantification/')

#organism = 'yeast'
organism = 'mouse'

fold = '10x'
prefix = paste0(organism, '_', fold)

target = 'trinity'

# read ambiguous site
ambiguous_ref_path = paste0(organism, '_data/ambiguous_site_', organism, '_cdna.tsv')
ambiguous_ref = read.table(ambiguous_ref_path, header=T, sep='\t', stringsAsFactors=F)

ambiguous_set_path = paste0(organism, '_data/ambiguous_site_', organism, '_', target, '_', fold,
                            '.tsv')
ambiguous_set = read.table(ambiguous_set_path, header=T, sep='\t',
                           stringsAsFactors=F)
rownames(ambiguous_ref) = ambiguous_ref$target_id
rownames(ambiguous_set) = ambiguous_set$target_id

# read ambiguous site (NR)
ambiguous_nr_path = paste0(organism, '_data/ambiguous_site_', organism, '_', target, '_', fold,
                           '_answer_label.tsv')
ambiguous_nr_set = read.table(ambiguous_nr_path, header=T, sep='\t', stringsAsFactors=F)
rownames(ambiguous_nr_set) = ambiguous_nr_set$target_id

# read assembly completeness
answer_label_path = paste0(organism, '_data/answer_label_', organism, '_', target, '_', fold, '.tsv')
answer_label = read.table(answer_label_path, header=F, stringsAsFactors=F, sep='\t') 
colnames(answer_label) = c('ref_id', 'set_id', 'orientation', 'ref_length', 'set_length',
                           'ref_cigar', 'set_cigar', 'ref_coverage', 'set_coverage')

answer_label$set_count = table(answer_label$set_id)[answer_label$set_id]
answer_label$category = 'others'
condition = answer_label$ref_coverage < 0.25 & answer_label$set_coverage < 0.25
answer_label[condition, 'category'] = 'noise'
condition = answer_label$ref_coverage >= 0.75 & answer_label$set_coverage >= 0.75
answer_label[condition, 'category'] = 'full-length'
condition = answer_label$ref_coverage < 0.25 & answer_label$set_coverage >= 0.75
answer_label[condition, 'category'] = 'fragmented'
condition = answer_label$ref_coverage >= 0.75 & answer_label$set_coverage < 0.25
answer_label[condition, 'category'] = 'over-extended'
condition = answer_label$ref_coverage >= 0.75 & answer_label$set_coverage >= 0.75
answer_label[condition, 'category'] = 'full-length'
condition = answer_label$set_count > 1
answer_label[condition, 'category'] = 'disaggregate'
answer_label$category = factor(answer_label$category)

# read disaggregate type
disaggregate_type_path = paste0(organism, '_data/disaggregate_type_', organism, '_', target, '_', fold,
                                '.tsv')
disaggregate_type = read.table(disaggregate_type_path, header=T, sep='\t', stringsAsFactors=F)
rownames(disaggregate_type) = disaggregate_type$target_id

# merge ambiguous sites
answer_label$ref_ambiguous_site = ambiguous_ref[answer_label$ref_id, 'ambiguous_sites']
answer_label$set_ambiguous_site = ambiguous_set[answer_label$set_id, 'ambiguous_sites']
answer_label$ref_ambiguous_cov = ambiguous_ref[answer_label$ref_id, 'align_coverage']
answer_label$set_ambiguous_cov = ambiguous_set[answer_label$set_id, 'align_coverage']

answer_label$ambiguous_type = 'unique_unique'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous'

####################################################################
# STOP HERE AND CHOOSE ONE OF A MERGING METHOD FOR AMBIGUOUS SITES #
####################################################################

# 1. merge ambiguous sites (with NR)
# 2. merge ambiguous sites (with NR, disaggregate)
# 3. merge ambiguous sites (with NR, disaggregate, disaggregate type-relative)
# 4. merge ambiguous sites (with NR, disaggregate, disaggregate type > 0, separate)
# 5. merge ambiguous sites (with NR, disaggregate, disaggregate type > 0, overlap)

###################################
# merge ambiguous sites (with NR) #
###################################
answer_label$set_nr_ambiguous_site = ambiguous_nr_set[answer_label$set_id, 'ambiguous_sites']
answer_label$ambiguous_type = 'unique_unique'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site > 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_ambiguous'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site == 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_unique'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site > 0 
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_ambiguous'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site == 0 
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_unique'

#################################################
# merge ambiguous sites (with NR, disaggregate) #
#################################################
answer_label$set_nr_ambiguous_site = ambiguous_nr_set[answer_label$set_id, 'ambiguous_sites']
answer_label$ambiguous_type = 'unique_unique_single'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site == 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'unique_unique_disaggregate'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category != 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_ambiguous_single'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_ambiguous_disaggregate'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category != 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_unique_single'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_unique_disaggregate'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0 &
            answer_label$category != 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique_single'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0 &
            answer_label$category == 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique_disaggregate'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category != 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_ambiguous_single'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_ambiguous_disaggregate'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category != 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_unique_single'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_unique_disaggregate'
#---------------------------------------------------------------------------------------------

#############################################################################
# merge ambiguous sites (with NR, disaggregate, disaggregate type-relative) #
#############################################################################
answer_label$separate = disaggregate_type[answer_label$set_id, 'separate']
answer_label$overlap = disaggregate_type[answer_label$set_id, 'overlap']
answer_label$set_nr_ambiguous_site = ambiguous_nr_set[answer_label$set_id, 'ambiguous_sites']
answer_label$ambiguous_type = 'unique_unique_single'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site == 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'unique_unique_disaggregate_separate'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site == 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate < answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'unique_unique_disaggregate_overlap'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site == 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate == answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'unique_unique_disaggregate_medium'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category != 'disaggregate' 
            answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_ambiguous_single'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_ambiguous_disaggregate_separate'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$separate < answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_ambiguous_disaggregate_overlap'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$separate == answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_ambiguous_disaggregate_medium'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category != 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_unique_single'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_unique_disaggregate_separate'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate < answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_unique_disaggregate_overlap'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate == answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_unique_disaggregate_medium'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0 &
            answer_label$category != 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique_single'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0 &
            answer_label$category == 'disaggregate' & answer_label$separate > answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique_disaggregate_separate'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0 &
            answer_label$category == 'disaggregate' & answer_label$separate < answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique_disaggregate_overlap'
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0 &
            answer_label$category == 'disaggregate' & answer_label$separate == answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique_disaggregate_medium'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category != 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_ambiguous_single'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_ambiguous_disaggregate_separate'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$separate < answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_ambiguous_disaggregate_overlap'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$separate == answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_ambiguous_disaggregate_medium'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category != 'disaggregate'
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_unique_single'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_unique_disaggregate_separate'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate < answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_unique_disaggregate_overlap'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate == answer_label$overlap
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_unique_disaggregate_medium'
#---------------------------------------------------------------------------------------------

##################################################################################
# merge ambiguous sites (with NR, disaggregate, disaggregate type > 0, separate) #
##################################################################################
answer_label$separate = disaggregate_type[answer_label$set_id, 'separate']
answer_label$overlap = disaggregate_type[answer_label$set_id, 'overlap']
answer_label$set_nr_ambiguous_site = ambiguous_nr_set[answer_label$set_id, 'ambiguous_sites']
answer_label$ambiguous_type = 'none_type'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site == 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > 0
answer_label[condition, 'ambiguous_type'] = 'unique_unique_disaggregate_separate_0'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_ambiguous_disaggregate_separate_0'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_unique_disaggregate_separate_0'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0 &
            answer_label$category == 'disaggregate' & answer_label$separate > 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique_disaggregate_separate_0'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > 0
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_ambiguous_disaggregate_separate_0'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$separate > 0
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_unique_disaggregate_separate_0'
#---------------------------------------------------------------------------------------------

################################################################################
# merge ambiguous sites (with NR, disaggregate, disaggregate type > 0, overlap #
################################################################################
answer_label$separate = disaggregate_type[answer_label$set_id, 'separate']
answer_label$overlap = disaggregate_type[answer_label$set_id, 'overlap']
answer_label$set_nr_ambiguous_site = ambiguous_nr_set[answer_label$set_id, 'ambiguous_sites']
answer_label$ambiguous_type = 'none_type'
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site == 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$overlap > 0
answer_label[condition, 'ambiguous_type'] = 'unique_unique_disaggregate_overlap_0'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$overlap > 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_ambiguous_disaggregate_overlap_0'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site > 0 & 
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$overlap > 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_ambiguous_unique_disaggregate_overlap_0'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site > 0 & answer_label$set_ambiguous_site == 0 &
            answer_label$category == 'disaggregate' & answer_label$overlap > 0
answer_label[condition, 'ambiguous_type'] = 'ambiguous_unique_disaggregate_overlap_0'
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site > 0 & answer_label$category == 'disaggregate' &
            answer_label$overlap > 0
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_ambiguous_disaggregate_overlap_0 '
#---------------------------------------------------------------------------------------------
condition = answer_label$ref_ambiguous_site == 0 & answer_label$set_ambiguous_site > 0 &
            answer_label$set_nr_ambiguous_site == 0 & answer_label$category == 'disaggregate' &
            answer_label$overlap > 0
answer_label[condition, 'ambiguous_type'] = 'unique_ambiguous_unique_disaggregate_overlap_0'
#---------------------------------------------------------------------------------------------

#####################################
# FINISHING MERGING AMBIGUOUS SITES #
#####################################

rownames(answer_label) = answer_label$ref_id
meta_data = answer_label

# some statistics
target_answer_label = subset(answer_label, answer_label$set_nr_ambiguous_site == 0)
print(table(target_answer_label$ambiguous_type))
target_answer_label = subset(answer_label, answer_label$set_nr_ambiguous_site != 0)
print(table(target_answer_label$ambiguous_type))

target_answer_label = subset(answer_label, answer_label$category == 'disaggregate')
print(table(target_answer_label$ambiguous_type))
target_answer_label = subset(answer_label, answer_label$category != 'disaggregate')
print(table(target_answer_label$ambiguous_type))

png(paste0('figures/', organism, '_ref_ambiguous_site.png'), width=500, height=400)
figure = ggplot(answer_label, aes(factor(category), ref_ambiguous_site))
figure + geom_boxplot(outlier.shape = NA) + ylim(0, 5) + xlab('category') + 
         ylab('reference ambiguous site') + ggtitle(paste0(organism, ' (reference)'))
dev.off()

png(paste0('figures/', organism, '_set_ambiguous_site.png'), width=500, height=400)
figure = ggplot(answer_label, aes(factor(category), set_ambiguous_site))
figure + geom_boxplot(outlier.shape = NA) + ylim(0, 5) + xlab('category') + 
    ylab('Trinity ambiguous site') + ggtitle(paste0(organism, ' (Trinity)'))
dev.off()

png(paste0('figures/', organism, '_category_count.png'), width=500, height=400)
figure = ggplot(answer_label, aes(x=category, fill=category))
figure + geom_bar()
dev.off()

