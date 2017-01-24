read_answer_label <- function(file_path){
    answer_label <- read.table(file_path, header=F, stringsAsFactors=F, sep='\t') 
    colnames(answer_label) <- c('ref_id', 'set_id', 'orientation', 'ref_length', 'set_length',
                                'ref_cigar', 'set_cigar', 'ref_coverage', 'set_coverage')
    set_count <- table(answer_label$set_id)[answer_label$set_id]
    answer_label$set_count <- as.vector(set_count)
    return(answer_label)
}

read_ambiguous_sites  <- function(file_path, answer_label, type){
    ambiguous = read.table(file_path, header=T, sep='\t', stringsAsFactors=F)
    rownames(ambiguous) = ambiguous$target_id
    if(type == 'ref' || type == 'reference'){
        answer_label$ref_ambiguous_site = ambiguous[answer_label$ref_id, 'ambiguous_sites']
        answer_label$ref_ambiguous_cov = ambiguous[answer_label$ref_id, 'align_coverage']
    } else if (type == 'set' || type == 'trinity'){
        answer_label$set_ambiguous_site = ambiguous[answer_label$set_id, 'ambiguous_sites']
        answer_label$set_ambiguous_cov = ambiguous[answer_label$set_id, 'align_coverage']
    } else if (type == 'nr' || type == 'set_nr' || type == 'trinity_answer_label') {
        answer_label$nr_ambiguous_site = ambiguous[answer_label$set_id, 'ambiguous_sites']
        answer_label$nr_ambiguous_cov = ambiguous[answer_label$set_id, 'align_coverage']
    }
    return(answer_label)
}

read_disaggregate_type <- function(file_path, answer_label){
    disaggregate_type = read.table(file_path, header=T, sep='\t', stringsAsFactors=F)
    rownames(disaggregate_type) = disaggregate_type$target_id
    answer_label$disaggregate_type_unique = disaggregate_type[answer_label$set_id, 'unique']
    answer_label$disaggregate_type_separate = disaggregate_type[answer_label$set_id, 'separate']
    answer_label$disaggregate_type_overlap = disaggregate_type[answer_label$set_id, 'overlap']    
    
    return(answer_label)
}

quantification_scatter_plot <- function(target_meta, figure_file, p_cor_1, s_cor_1, p_cor_2, s_cor_2){
    figure_1 = ggplot(target_meta, aes(ans_1, set_count_1, shape='a', colour="#FF9999"))
    figure_1 = figure_1 + geom_point() + 
               annotate("text", x=0, y=Inf, hjust = 0, vjust=1.25, size=3.5,
                        label=paste0('transcript count: ', nrow(target_meta), '\n', "Pearson's R: ", p_cor_1, '\n', "Spearman's R: ", s_cor_1, sep='')) +
               scale_shape_discrete(solid=F, guide=F) + 
               scale_colour_discrete(guide=F) + 
               xlab('exact count') +
               ylab('estimate count') +
               xlim(0, max(quantile(target_meta$ans_1, 0.975),  quantile(target_meta$set_count_1, 0.975))) + 
               ylim(0, max(quantile(target_meta$ans_1, 0.975),  quantile(target_meta$set_count_1, 0.975))) + 
               geom_smooth(method="lm",se=FALSE, color='black')
    
    figure_2 = ggplot(target_meta, aes(ans_2, set_count_2, shape='a', colour="#FF9999"))
    figure_2 = figure_2 + geom_point() + 
               annotate("text", x=0, y=Inf, hjust = 0, vjust=1.25, size=3.5,
                        label=paste0('transcript count: ', nrow(target_meta), '\n', "Pearson's R: ", p_cor_2, '\n', "Spearman's R: ", s_cor_2, sep='')) +
               scale_shape_discrete(solid=F, guide=F) + 
               scale_colour_discrete(guide=F) + 
               xlab('exact count') +
               ylab('estimate count') +
               xlim(0, max(quantile(target_meta$ans_2, 0.975),  quantile(target_meta$set_count_2, 0.975))) + 
               ylim(0, max(quantile(target_meta$ans_2, 0.975),  quantile(target_meta$set_count_2, 0.975))) + 
               geom_smooth(method="lm",se=FALSE, color='black')

    jpeg(figure_file, res=350, width=3000, height=1500)
    multiplot(figure_1, figure_2, cols=2)
    dev.off()
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
    grid.newpage()
    grid.draw(combined)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}
