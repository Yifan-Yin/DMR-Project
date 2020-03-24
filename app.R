## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------------------
library(tidyverse)
library(egg)
library(here)


## --------------------------------------------------------------------------
# sample data
pDat <- readRDS(here('data', 'pDat.rds'))
# methylation data
betas <- readRDS(here('data', 'betas.rds'))
# annotation
anno <- readRDS(here('data', 'annotation.rds'))
#color code
color_code <- readRDS(here::here('data', '2_3_color_code.rds'))
color_code_tissue <- setNames(color_code$Colors_Tissue, color_code$label)
colors <- color_code_tissue[unique(pDat$Tissue)]

#  transcript models for exapnded view
anno_annotatr <- readRDS(here('data', 'annotation_annotatr.rds'))


## --------------------------------------------------------------------------
# rename column names
rename_columns <- function(x) {
  dplyr::rename(
    x,
    
    # columns to rename
    CpG = cpg,
    Chr = chr,
    Start = start,
    End = end,
    `Relation to CpG Islands` = cpg_id,
    `Enhancers (FANTOM5)` = enhancers_id,
    `Gene Symbol` = genes_symbol,
    `Relation to Transcript Features` = genes_id,
    `Relation to Transcripts` = genes_tx_id,
    `Relation to Placental PMD Regions` = pmd_id,
    `Relation to Imprinted DMRs (general)` = imprinted_dmr_general,
    `Relation to Imprinted DMRs (placenta)` = imprinted_dmr_placenta
  )
}


## --------------------------------------------------------------------------
#column names 
column_names <- anno %>% 
  select(cpg, chr, start, end, 
         cpg_id, enhancers_id, genes_symbol, genes_tx_id, genes_id, pmd_id, 
         imprinted_dmr_general, imprinted_dmr_placenta) %>%
  rename_columns() %>%
  names()


## --------------------------------------------------------------------------
#list of genes for search list
gene <- tibble(gene = str_split(anno$genes_symbol, ', ') %>%
                 unlist() %>%
                 unique() %>%
                 na.omit() %>%
                 sort())

#pull out cpg name for search list use
cpg <- unique(anno$cpg) 
cpg <- cpg[cpg %in% rownames(betas)] %>% as_tibble()


## --------------------------------------------------------------------------
#produce a datatable which shows cpg annotation
make_datatable <- function(cpg_name = NULL, gene_name = NULL, 
                           Chr = NULL, Start = NULL, End = NULL){
  if (!is.null(cpg_name)) { #get cpg annotation 
    anno_cpg <- filter(anno, cpg %in% cpg_name) %>% arrange(desc(start))
  }
  else if (!is.null(gene_name)) {
    #get cpg sites related to input genes
    cpg_name <- find_cpg_ids(gene_name)
    
    #get cpg annotation 
    anno_cpg <- filter(anno, cpg %in% cpg_name)%>% arrange(desc(start))
  }
  
  else if (!is.null(Chr) & !is.null(Start) & !is.null(End))
  {
    #get cpg sites within a given region
    cpg_name <- find_cpg_from_region(Chr, Start, End)
    
    #get cpg annotation 
    anno_cpg <- filter(anno, cpg %in% cpg_name)%>% arrange(desc(start)) 
  }
  
  anno_cpg <- anno_cpg %>% 
    select(cpg, chr, start, end, 
           cpg_id, enhancers_id, 
           genes_symbol, genes_tx_id, genes_id, 
           pmd_id,
           imprinted_dmr_general, imprinted_dmr_placenta) %>%
    rename_columns()
}



## --------------------------------------------------------------------------
annoPlot_with_tracks <- function(cpg_name = NULL, 
                                 gene_name = NULL, 
                                 Chr = NULL, Start = NULL, 
                                 End = NULL, cpg_number = 500, 
                                 first_cpg = NULL, end_cpg = NULL, 
                                 Tissue_type = unique(pDat$Tissue),
                                 Sex_type = c('F', 'M'),
                                 Trimester_type = c('Third', 'First'),
                                 plot_heights = c(36, 9, 6, 5)){
  
  ######## FILTER ANNOTATION TO RELEVANT CPGS, DEPENDING ON INPUT #################
  if (!is.null(gene_name) & is.null(first_cpg) & is.null(end_cpg)){
    
    # If a gene symbol is supplied,
    # filter to gene in annotation
    anno_gene <- anno %>%
      filter(grepl(paste0('(?<!-)\\b', gene_name, '\\b(?!-)'), genes_symbol, perl = TRUE),
             cpg %in% rownames(betas)) %>%
      arrange(desc(start)) %>% slice(1:cpg_number)
    
    # stop if gene symbol is invalid
    if (nrow(anno_gene) == 0) {
      stop('No cpgs found, maybe they were filtered out.')
    }
  } else if (!is.null(gene_name) & !is.null(first_cpg) & !is.null(end_cpg)){
    
    # If a gene symbol is supplied,
    # filter to gene in annotation
    anno_gene <- anno %>%
      filter(grepl(paste0('(?<!-)\\b', gene_name, '\\b(?!-)'), genes_symbol, perl = TRUE),
             cpg %in% rownames(betas)) %>%
      arrange(desc(start)) %>% slice(first_cpg: end_cpg)
    
    # stop if gene symbol is invalid
    if (nrow(anno_gene) == 0) {
      stop('No cpgs found, maybe they were filtered out.')
    }
  }  else if (!is.null(cpg_name)){
    
    # If a cpg is provided,
    # filter annotation to cpgs
    anno_gene <- anno %>%
      filter(cpg %in% cpg_name,
             cpg %in% rownames(betas)) %>% 
      arrange(chr, desc(start))
    
    # stop if no cpgs were found
    if (nrow(anno_gene) == 0) {
      stop('None of the selected CpGs exist in our processed data.')
    }
    
    if (length(cpg_name) < length(intersect(cpg_name, rownames(betas)))) {
      print(paste0('The following selected CpGs do not exist in our processed data:\n', 
                   setdiff(cpg_name, rownames(betas))))
    }
    
  } else if (!is.null(Chr) & !is.null(Start) & !is.null(End) & is.null(first_cpg) & is.null(end_cpg)){
    
    #get cpg sites within a given region
    cpg_name <- find_cpg_from_region(Chr, Start, End)
    # filter annotation to cpgs
    anno_gene <- anno %>%
      filter(cpg %in% cpg_name,
             cpg %in% rownames(betas)) %>% 
      arrange(chr, desc(start)) %>% slice(1:cpg_number)
    # stop if no cpgs were found
    if (nrow(anno_gene) == 0) {
      stop('No cpgs found in the selected region')
    }
  } else if (!is.null(Chr) & !is.null(Start) & !is.null(End) & !is.null(first_cpg) & !is.null(end_cpg)){
    
    #get cpg sites within a given region
    cpg_name <- find_cpg_from_region(Chr, Start, End)
    # filter annotation to cpgs
    anno_gene <- anno %>%
      filter(cpg %in% cpg_name,
             cpg %in% rownames(betas)) %>% 
      arrange(chr, desc(start)) %>% slice(first_cpg : end_cpg)
    
    # stop if no cpgs were found
    if (nrow(anno_gene) == 0) {
      stop('No cpgs found in the selected region')
    }
  } else {
    stop('Enter valid cpg(s) OR one specific gene.')
    
  }
  
  ######## SETUP INPUT TO PLOTS ###############
  #
  # Title
  #
  title <- paste0(
    unique(anno_gene$chr), ':',
    scales::comma(min(anno_gene$start)), '-',
    scales::comma(max(anno_gene$start))
  )
  
  if (is.null(cpg_name)){
    title <- paste0(gene_name, ', ', title)
  } 
  
  
  #
  # now we wrangle the betas/annotations/pdata information into a format usable for plotting
  #
  betas_gene <- t(betas[anno_gene$cpg,,drop = FALSE]*100) %>% 
    as_tibble %>% 
    mutate(deidentified_id = colnames(betas)) %>%
    
    # reshape into longer format
    pivot_longer(cols = -deidentified_id,
                 names_to = 'cpg', 
                 values_to = 'beta')  %>%
    
    # add tissue and trimester info
    left_join(pDat %>% select(deidentified_id, Trimester, Tissue, Sex), by = 'deidentified_id') %>%
    
    filter(Sex %in% Sex_type,
           Trimester %in% Trimester_type) %>%
    
    # calculate mean and sd for each cpg for each group
    group_by(Tissue, Trimester, cpg) %>%
    summarize(mean = mean(beta),
              sd = sd(beta)) %>%
    mutate(lower = mean-sd, upper = mean+sd) %>%
    
    # add cpg info
    left_join(anno_gene, by = 'cpg') %>%
    
    # order on cpg position
    ungroup() %>%
    arrange(desc(start)) %>%
    mutate(cpg = factor(cpg, levels = unique(cpg)),
           Trimester_Tissue = paste0(Trimester, ' - ', Tissue)) %>%
    filter(Tissue %in% Tissue_type)
  
  ##### GENERATE INDIVIDUAL PLOT TRACKS #####
  font_size <- 14
  point_size <- 3
  line_size <- 1.25
  
  # Generate  methylation track
  p1 <- betas_gene %>%
    ggplot() +
    
    geom_blank(aes(x = cpg)) +
    # geom_ribbon is to generate the alternating shaded background
    #geom_ribbon(aes(x = as.numeric(cpg), ymin = 0, ymax = 25),
    #            fill = 'grey', alpha = 0.15)+
    #geom_ribbon(aes(x = as.numeric(cpg), ymin = 50, ymax = 75),
    #            fill = 'grey', alpha = 0.15)+
    
    # geom_linerange and geom_point is for the actual methylation data
    geom_linerange(alpha = 0.5, size = line_size, 
                   aes(x = cpg, ymin =lower, ymax = upper, color = Tissue),
                   show.legend = FALSE) +
    geom_point(alpha = 0.75, aes(x = cpg, y = mean, color = Tissue), size = point_size) +
    
    # we want to facet by trimester
    facet_wrap(vars(Trimester), ncol = 1, labeller = label_both, scales = 'free') +
    
    # cosmetics
    theme_bw(base_size = font_size) +
    theme(axis.text.x = element_blank(), 
          #axis.text.x = element_text(angle = 90),
          axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0, size = 16),
          
          legend.text = element_text(size = 14),
          
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          
          panel.spacing.y = unit(0.1, 'cm'),
          
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold', hjust = 0),
          
          plot.margin = unit(c(0, 2, 0, 2), 'pt'),) +
    scale_y_continuous(limits = c(0, 100), 
                       expand = c(0, 0), 
                       breaks = c(0, 50, 100),
                       labels = function(x)paste0(x, '%')) +
    scale_x_discrete(expand = c(0,0.5))+
    scale_color_manual(values= color_code_tissue[unique(betas_gene$Tissue)],
                       guide = guide_legend(override.aes = list(size = 4)),
                       labels = function(x)gsub(' cs', '', x)) +
    labs(y = 'DNA\nmethylation', 
         x = '', 
         
         color = '', 
         title = title)
  
  
  p_dmc <- betas_gene %>%
    
    # tidy data
    select(cpg, `First_Endothelial cs`:`Third_Trophoblasts cs`) %>%
    distinct() %>%
    pivot_longer(cols = -cpg,
                 names_to = 'Group',
                 values_to = 'Significant') %>%
    mutate(facet_title = 'Differential methylation by cell type') %>%
    separate(Group, c('Trimester', 'Tissue'), '_') %>%
    
    # convert axes to numeric for geom_path to work
    mutate(cpg_num = as.numeric(cpg),
           Tissue_num = as.numeric(as.factor(Tissue)),
           
           # Add tissue-specific color when significant
           Significant_num = if_else(Significant, Tissue_num, NA_real_)) %>% 
    filter(Trimester %in% Trimester_type) %>%
    
    {
      ggplot(data = ., aes(y = Significant_num, color = Tissue)) +
        geom_blank(aes(x = cpg)) +
        geom_point(data = . %>% 
                     filter(!is.na(Significant_num)), size = point_size,
                   aes(x =cpg_num)) +
        geom_path(na.rm = TRUE, size = line_size-0.25, alpha = 1, aes(x =cpg_num)) +
        facet_grid(rows = vars(Trimester),
                   cols = vars(facet_title),
                   labeller =  labeller(Trimester = function(x)paste0('Trimester:\n', x)),
                   switch = 'y') +
        theme_bw(base_size = font_size) +
        theme(axis.text.x = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect(fill = '#f7f7f7'),
              
              #axis.text.x = element_text(angle = 90, vjust = 0),
              axis.ticks = element_blank(),
              axis.text.y = element_text(vjust = 0.5),
              
              legend.text = element_text(size = 14),
              
              plot.margin = unit(c(0, 2, 0, 2), 'pt'),
              plot.background = element_rect(fill = '#f7f7f7'),
              
              strip.background = element_blank(),
              strip.text.y = element_text(face = 'bold', vjust = 0.5, hjust = 0.5, angle = 180),
              strip.text.x = element_text(face = 'bold', hjust = 0),
              strip.placement = 'outside') +
        scale_color_manual(values = color_code_tissue[unique(pDat$Tissue)],
                           na.value = '#f7f7f7',
                           guide = 'none') +
        scale_y_continuous(expand = c(0.15,0.15),
                           breaks = c(1, 2, 3, 4, 5), 
                           labels = pDat$Tissue %>% 
                             as.factor() %>% 
                             levels() %>%
                             gsub(' cs', '', .)) +
        scale_x_discrete(expand = c(0,0.5)) +
        labs(x = '', y = '')
      
    }
  
  # plot cpg annotations
  p2 <- betas_gene %>% 
    
    # First we need to process the annotation data
    mutate(# absent/presence for different genomic elements
      enhancer = !is.na(enhancers_id),
      pmd = !is.na(pmd_id)) %>%
    mutate_if(is.logical, as.character) %>%
    select(cpg, enhancer, pmd, 
           imprinted_dmr_general, imprinted_dmr_placenta, cpg_id, start) %>%
    
    # reshape
    pivot_longer(cols = c(-cpg, -start),
                 names_to = 'cpg_element', 
                 values_to = 'presence') %>%
    distinct() %>%
    
    # arrange plot orde for each track
    arrange(desc(start)) %>%
    mutate(cpg = factor(cpg, levels = unique(cpg)),
           cpg_element = factor(cpg_element, 
                                levels = c('enhancer', 'pmd', 'imprinted_dmr_general',
                                           'imprinted_dmr_placenta', 'cpg_id')),
           presence = factor(ifelse(presence == "TRUE", 'Present',
                                    ifelse(presence == "FALSE", 'Absent', presence)),
                             levels = c('sea', 'shelf', 'shore', 'island',
                                        'Present', 'Absent', ' ', '  ')),
           # for facet label
           group = 'Annotations') %>%
    
    # plot code
    ggplot(aes(x = cpg, y = cpg_element, fill = presence)) +
    geom_tile(color = 'white') +
    guides(fill=guide_legend(ncol=2, override.aes = list(colour = 'white'))) +
    facet_wrap(vars(group), ncol = 2) +
    theme_bw(base_size = font_size) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(vjust = 0.5),
          
          legend.text = element_text(size = 14),
          
          plot.margin = unit(c(0, 2, 0, 2), 'pt'),
          
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold', hjust = 0)) +
    scale_fill_manual(values = c('Present' = '#cccccc', 'Absent' = '#f7f7f7', 
                                 ' ' = 'white', '  ' = 'white',
                                 'sea' = '#ffffcc', 'shelf' = '#a1dab4', 'shore' = '#41b6c4',
                                 'island' = '#225ea8'), drop = F,
                      labels = stringr::str_to_title) +
    scale_y_discrete(expand = c(0,0),
                     labels = function(x)gsub('cpg_id', 'Relation to CpG Islands',
                                              gsub('imprinted_dmr_placenta', 'Imprinted DMR (Placenta)',
                                                   gsub('imprinted_dmr_general', 'Imprinted DMR (General)',
                                                        gsub('pmd', 'Placental PMD',
                                                             gsub('enhancer', 'Enhancer', 
                                                                  x)))))) +
    scale_x_discrete(expand = c(0,0.2)) +
    labs(x = '', y = '', fill = '')
  
  # transcript annotations
  
  # cpgs can have more than one mapping to a gene element per transcript
  # cpgs mapped to intronexonboundaries can also be part of exons and introns
  # 5 UTR and 3 UTR can also overlap with exons and introns
  
  # If there are multiple mappings then only the following will be displayed in order of priority:
  # 5UTR / 3UTR, exon / intron, intronexonboundary
  
  p3 <- betas_gene %>% 
    select(cpg, end, genes_id, genes_tx_id, genes_symbol) %>%
    distinct() %>%
    
    # split comma separate list and put each entry into a new row
    mutate(genes_tx_id = str_split(genes_tx_id, ', '),
           genes_id = str_split(genes_id, ', '),
           genes_symbol = str_split(genes_symbol, ', ')) %>%
    unnest(c(genes_tx_id, genes_id, genes_symbol)) %>%
    
    # paste together transcript ID and genes symbol
    mutate(genes_tx_id_symbol = if_else(genes_symbol != 'no_associated_gene', 
                                        paste0(genes_tx_id, ' (', genes_symbol, ')'),
                                        genes_tx_id)) %>%
    
    ### ORDER TRANSCRIPTS IN PLOT 
    #
    # note that the first level appears on the bottom of the graph (y = 0 if it was continuous)
    
    # calculate earliest end of each symbol
    group_by(genes_symbol) %>%
    mutate(genes_symbol_end = min(end)) %>%
    ungroup() %>%
    
    # reorder gene symbols by their earliest end position
    mutate(genes_symbol = fct_reorder(genes_symbol, genes_symbol_end)) %>%
    
    # reorder transcripts based on gene symbol levels and then by their end position
    arrange(genes_symbol, end) %>%
    mutate(genes_tx_id_symbol = 
             factor(genes_tx_id_symbol, levels = unique(genes_tx_id_symbol)) %>%
             fct_explicit_na()) %>%
    
    ####
    # fill all combinations of cpg and gene_tx_id with NA if missing, so that they show in plot
    complete(cpg, genes_tx_id_symbol) %>%
    
    # filter out to unique mappings using the priority described above
    # strategy is to order by factor level and take the first occurence
    mutate(genes_id = factor(genes_id, 
                             levels = c('1to5kb', 'promoter', '5UTR', '3UTR', 'exon',
                                        'intron', 'intronexonboundary', 'intergenic'))) %>%
    
    group_by(cpg, genes_tx_id_symbol) %>%
    arrange(genes_id)  %>%
    dplyr::slice(1) %>%
    
    # reorder factor levels for appearance in plot legend
    mutate(genes_id = factor(genes_id, 
                             levels = c('1to5kb', 'promoter', '5UTR', 'exon', 'intron',
                                        'intronexonboundary', '3UTR', 'intergenic')),
           group = 'Transcripts') %>%
    
    #plot
    ggplot(aes(x = cpg, y = genes_tx_id_symbol, fill = genes_id)) +
    geom_tile(color = 'white', alpha = 0.6) +
    facet_wrap(vars(group), ncol = 1) +
    guides(fill=guide_legend(ncol=1, override.aes = list(alpha=0.7, color = 'white'))) +
    theme_bw(base_size = font_size) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          axis.text.y = element_text(vjust = 0.5),
          axis.ticks = element_blank(),
          
          plot.margin = unit(c(0, 2, 0, 2), 'pt'),
          
          legend.text = element_text(size = 14),
          
          panel.grid = element_blank(),
          
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold', hjust = 0)) +
    scale_fill_brewer(na.value = '#f7f7f7', breaks = c('1to5kb', 'promoter', '5UTR', 'exon', 
                                                       'intron', 'intronexonboundary', '3UTR',
                                                       'intergenic'),
                      palette = 'Paired', direction = -1) +
    scale_y_discrete(expand = c(0,0), 
                     labels = function(x)(if_else(x == '(Missing)', '', x))) +
    scale_x_discrete(expand = c(0,0.2)) +
    coord_cartesian(clip = "off") +
    labs(x = '', y = '', fill = '')
  
  # remove cpg names iftoo many cpgs are displayed
  if (nrow(anno_gene > 50)) {
    p3 <- p3 + theme(axis.text.x = element_blank())
    
  }
  # the size of each plot depends on the number of unique transcripts
  # More transcripts means a longer transcript track
  # so the plot height depends on the number of transcripts
  
  # number of transcripts
  n_trans <- betas_gene$genes_tx_id %>% str_split(', ') %>% unlist() %>% unique %>% length()
  
  
  # maximum size of transcript plot is "25"
  n_trans <- min(c(n_trans, 25))
  
  # size of p1
  p1_h <- 40
  
  # combine tracks plot
  if (is.null(plot_heights)){
    ggarrange(p1, p_dmc, p2, p3, ncol = 1, heights = c(p1_h-(n_trans-6), 9, 6, n_trans))
  } else {
    ggarrange(p1, p_dmc, p2, p3, ncol = 1, heights = plot_heights)
  }
  
  
}



## --------------------------------------------------------------------------
find_cpg_ids <- function(gene_name){
  cpg_name <- anno %>%
    filter(cpg %in% rownames(betas),
           grepl(paste0('(?<!-)\\b', gene_name, '\\b(?!-)'), genes_symbol, perl = TRUE)) %>%
    pull(cpg) %>% 
    unique()
  
  cpg_name
}


## --------------------------------------------------------------------------
find_cpg_from_region <- function(Chr, Start, End){
  
  cpg_name <- anno %>% 
    filter(chr == Chr, start > Start, end < End) %>%
    pull(cpg) %>%
    unique()
  
  cpg_name <- cpg_name[cpg_name %in% rownames(betas)]
  
  cpg_name
}



## --------------------------------------------------------------------------
sample_info_cpg <- function(cpg_name){
  
  cpg <- betas[cpg_name,]
  
  #when have only one cpg input
  if (length(cpg_name) ==1){
    cpg <- as.data.frame(cpg) %>% t()
    cpg <- cpg[,colnames(cpg) %in% pDat$deidentified_id]
    cpg <- as.data.frame(cpg)
    colnames(cpg) <- cpg_name
    cpg <- cpg %>% cbind(pDat) %>% gather(cpg, betas, contains('cg'))
  }
  
  #when have more than one cpg inputs
  else{
    cpg <- cpg[,colnames(cpg) %in% pDat$deidentified_id]
    cpg <- cpg %>% t() %>% cbind(pDat) %>% gather(cpg, betas, contains('cg'))
  }
  cpg <- cpg %>% 
    left_join(anno, by = 'cpg') %>% 
    arrange(desc(start)) %>% 
    mutate(cpg = factor(cpg, levels = unique(cpg)))
  cpg
}


## --------------------------------------------------------------------------
#boxplots which show beta values based on cell types
boxplot_selected<- function(cpg_name, sex_type = c('F', 'M')){
  
  #given cpg, get sample betas, cell type  and trimester
  cpg <- sample_info_cpg(cpg_name) %>%
    filter(Sex %in% sex_type)
  
  ggplot(cpg, aes(x=Trimester, y = betas, color = Tissue)) +
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(alpha = 0.9, show.legend = F, 
                position = position_jitterdodge(jitter.width = 0.05)) +
    facet_wrap(.~cpg)+ 
    theme_bw(base_size = 16)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0))+
    labs(x = ' ', y ='DNA\nmethylation', color = '') +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) 
}


## --------------------------------------------------------------------------
expanded_view <- function(gene_name = NULL, 
                          
                          start_extend = 2500, # how far to extend the plot left, in genomic coordinates 
                          end_extend = 2500,   # how far to extend right
                          
                          show_transcript_gene = TRUE, # show transcripts associated gene symbol
                          breaks_width = 15000, # determines the spacing of x axis breaks
                          
                          highlight_region = NULL, # a df with start end columns, indicating regions to highlight
                          
                          trimester = c('Third', 'First'), # can be Third, First, or Both
                          sex = c('F', 'M'), # can be Male, Female, or Both
                          tissue = c("Endothelial cs", "Villi", "Hofbauer cs", 
                                     "Trophoblasts cs", "Stromal cs"),
                          plot_sizes = c(3,1) # a 2-element numeric vector indicating the relative plot heights 
                          # of the methylation plot to the transcript plot
                          
                          # required components:
                          # pDat with columns Sample_Name, Trimester, Sex, and Tissue
                          # betas_filt data.frame with column names == pDat$Sample_Name
                          # anno_annotatr data frame that is contains a column cpg, genes_symbol, chr, start, end
                          # the cpgs in anno_annotatr and betas_filt need to be exactly the same, and same order
                          # annotation data.frame that I built from annotatr, and slightly cleaned
                          
) {
  
  # subset to relevant samples
  pDat <- pDat %>%
    filter(Trimester %in% trimester,
           Sex %in% sex,
           Tissue %in% tissue)
  
  # subset to cpgs associated with gene
  anno_gene_subset <- anno %>%
    filter(grepl(paste0('\\<', gene_name, '\\>'), genes_symbol)) %>% 
    dplyr::select(chr, start, end) 
  
  c <- unique(anno_gene_subset$chr)
  s <- min(anno_gene_subset$start) # first cpg
  e <- max(anno_gene_subset$start) # last cpg
  s_extend <- s - start_extend
  e_extend <- e + end_extend
  
  region_of_interest <- anno %>%
    filter(cpg %in% rownames(betas),
           chr == c, start > s_extend, end < e_extend)
  
  # subset betas
  b <- betas[region_of_interest$cpg, pDat$deidentified_id] %>%
    as.data.frame() %>%
    bind_cols(cpg = rownames(.), .) %>%
    pivot_longer(cols = -cpg,
                 names_to = 'deidentified_id',
                 values_to = 'beta') %>%
    left_join(anno, by = 'cpg') %>% 
    right_join(pDat, by = 'deidentified_id') 
  
  
  p1 <- b %>%
    
    # calculate mean beta for each cpg for each tissue
    group_by(cpg, Tissue) %>%
    summarize(mean = mean(beta), start = dplyr::first(start)) %>%
    
    # plot methylation means
    ggplot(aes(x = start)) +
    geom_linerange(ymin = 0, aes(ymax = mean, color = Tissue)) +
    
    facet_wrap(vars(Tissue), ncol = 1, strip.position = 'left',
               labeller = labeller(Tissue = function(x)gsub(' cs', '', x))) +
    
    theme_bw(base_size = 16) +
    theme(strip.background = element_blank(),
          strip.text.y = element_text(angle = 180, hjust = 1),
          strip.placement = 'outside',
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_line(size = 0.5),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0, 'lines'),
          #panel.grid = element_blank()
    ) +
    scale_y_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0,0), labels = NULL) +
    scale_x_continuous(limits = c(s_extend, e_extend), expand = c(0,0), 
                       breaks = scales::breaks_width(15000),
                       labels = scales::number) +
    scale_color_manual(values = colors, guide = 'none') +
    labs(y = '', x = '', title = gene_name, 
         subtitle = paste0(c, ':', scales::comma(s), '-', scales::comma(e)));p1
  
  p2 <- anno_annotatr %>%
    filter(seqnames == as.character(c), start > s_extend, end < e_extend,
           type %in% c('exons', 'introns', 'promoters')) %>%
    
    mutate(symbol = ifelse(is.na(symbol), '', symbol)) %>%
    ggplot(aes(x = start, xend = end)) +
    geom_segment(y= 0.5, yend = 0.5, aes(size = type, color = type))  +
    facet_wrap(~symbol + tx_id, ncol = 1, strip.position = 'left') +
    theme_bw(base_size = 16)+
    theme(strip.background = element_blank(),
          strip.text.y  = element_text(angle = 180),
          strip.placement = 'outside',
          legend.position = 'bottom',
          plot.title =  element_blank(),
          plot.subtitle = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(size = 0.5),
          panel.spacing = unit(0, 'lines'),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.border = element_blank(),
          legend.key.size = unit(0.5, 'cm'),
          legend.spacing = unit(0.25, 'cm')) +
    scale_size_manual(values = c('exons' = 7, 'introns' = 1, 'intronexonboundaries' = 3,
                                 '5UTRs' = 5, '3UTRs' = 6,
                                 '1to5kb' = 2,
                                 'promoters' = 2,
                                 
                                 'cpg_inter' = 1,
                                 'cpg_islands' = 5,
                                 'cpg_shores' = 4,
                                 'cpg_shelves' = 3)) +
    scale_color_manual(values = c('exons' = 'black', 'introns' = 'black', 'intronexonboundaries' = 'black',
                                  '5UTRs' = 'grey', '3UTRs' = 'grey',
                                  '1to5kb' = 'grey',
                                  'promoters' = 'forestgreen',
                                  
                                  'cpg_inter' = 1,
                                  'cpg_islands' = 5,
                                  'cpg_shores' = 4,
                                  'cpg_shelves' = 3)) +
    scale_x_continuous(limits = c(s_extend,e_extend), 
                       expand = c(0,0), 
                       labels = scales::number,
                       breaks = scales::breaks_width(15000)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = 'position', size = '', color = '');p2
  
  if (!is.null(highlight_region)) {
    p1 <- p1 + 
      geom_rect(data = highlight_region,
                aes(xmin = start, xmax = end),
                ymin = 0, ymax = 1,
                alpha = 0.1)
    p2 <- p2 + 
      geom_rect(data = highlight_region,
                aes(xmin = start, xmax = end),
                ymin = 0, ymax = 1,
                alpha = 0.1)
  }
  egg::ggarrange(p1,p2, ncol = 1, heights = plot_sizes) 
  
}

## --------------------------------------------------------------------------
#install.packages('rsconnect')
#install.packages('shinyjs')
library(rsconnect)
library(shiny)
rsconnect::setAccountInfo(name='yifan7', token='52E1037550A7AA68252579A98EE8AA38', secret='+NlrJysSVZHa5LQQb9slmYRDXExPfjF0QbxMjUH1')
library(DT)
library(shinyjs)
library(shinythemes)


## --------------------------------------------------------------------------
ui <- fluidPage(
  theme = shinythemes::shinytheme("lumen"),
  
  useShinyjs(),
  #main title
  titlePanel("Welcome to the Placental Cell Methylome Browser!"),
  p('Scroll down to see detailed instructions.'),
  
  #UI layout
  sidebarLayout(
    sidebarPanel(
      
      # relative width out of 12 units
      width = 3,
      
      tags$h4(
        tags$strong('Plot input:')),
      
      p('Select a gene, specific CpGs, or a genomic region:', 
        span(tags$a(href = '#footer_1', tags$sup(1), style = 'color:blue'))),
      
      tabsetPanel(type = 'tabs', id = 'tabs',
                  
                  #page of gene input 
                  tabPanel('Gene Input',
                           
                           #gene search list
                           selectizeInput(
                             'gene_name', '',
                             choices = NULL,
                             multiple = FALSE,
                             options = list(placeholder='Enter a gene name',
                                            onInitialize = I('function() { this.setValue(""); }'))),
                           #action button to submit and reset input
                           actionButton("submit_gene_input", "Plot"),
                           actionButton('reset_gene_input','Clear')),
                  
                  #page of cpgs input
                  #tabPanel('CpG Input',
                  
                  #cpg search list
                  #         selectizeInput('cpg_name', '',
                  #                        choices = NULL,
                  #                        multiple = TRUE,
                  #                        options = list(placeholder='Enter CpG site names' )),
                  
                  #action button to submit and reset input
                  #         actionButton("submit_cpg_input", "Plot"),
                  #         actionButton('reset_cpg_input','Clear')),
                  
                  #page of region input
                  tabPanel('Region Input',
                           
                           #chromosome 
                           selectizeInput('chr','Chromosome',
                                          choices = NULL,
                                          multiple = FALSE,
                                          options = list(placeholder = 'Enter a chromosome')),
                           
                           #start point
                           numericInput('start', 'Start', value = NULL),
                           #end point
                           numericInput('end', 'End', value = NULL),
                           
                           #action button to submit and reset input
                           actionButton("submit_region_input", "Plot"),
                           actionButton('reset_region_input','Clear'))),
      hr(),
      
      tags$h4(
        tags$strong('Customize plot')),
      checkboxGroupInput('checkbox_tissue',
                         'Show methylation for:',
                         setNames(unique(pDat$Tissue), gsub(' cs', '', unique(pDat$Tissue))),
                         selected = unique(pDat$Tissue)),
      
      
      checkboxGroupInput('checkbox_sex',
                         'Include:',
                         setNames(c('F', 'M'), c('Female', 'Male')),
                         selected = unique(pDat$Sex)),
      
      checkboxGroupInput('checkbox_trimester',
                         'Include:',
                         setNames(c('First', 'Third'), c('First Trimester', 'Third Trimester (Term)')),
                         selected = unique(pDat$Trimester)),
      
      sliderInput('range_cpg', 
                  'Select number of CpGs to plot (only applies to Condensed View):',
                  min = 1, max = 500, value = c(500)),
      actionButton('submit_cpg_set','Submit'),
      br(),
      hr(),
      
      tags$h4(
        tags$strong('Customize CpG annotation table')),
      
      checkboxGroupInput('checkbox',
                         'Show columns:',
                         
                         choices = names(rename_columns(anno))[-13:-20],
                         inline = FALSE,
                         selected = names(rename_columns(anno))[-13:-20])),
    
    mainPanel(
      tabsetPanel(
        type = 'tabs',
        
        # First tab
        tabPanel("Condensed View",
                 
                 fluidRow(
                   br(),
                   tags$h4(
                     tags$strong(
                       'Adjust heights of individual data tracks:'
                     )
                   ),
                   column(
                     12,
                     wellPanel(
                       fluidRow(
                         column(
                           3, 
                           sliderInput(
                             'condensed_plot_height_1',
                             'Methylation',
                             min = 0, max = 50, value = c(35), step = 1)
                         ),
                         column(
                           3,
                           sliderInput(
                             'condensed_plot_height_2',
                             'Differential Methylation',
                             min = 0, max = 50, value = c(9), step = 1)),
                         column(
                           3, 
                           sliderInput(
                             'condensed_plot_height_3',
                             'Genomic features',
                             min = 0, max = 50, value = c(6), step = 1)),
                         column(
                           3, 
                           sliderInput(
                             'condensed_plot_height_4',
                             'Transcripts',
                             min = 0, max = 50, value = c(5), step = 1)
                         )
                       )
                     ))
                 ),
                 
                 br(),
                 textOutput('text'),
                 plotOutput('annoplot', height = 1100)),
        
        tabPanel("Expanded View",
                 
                 fluidRow(
                   br(),
                   tags$h4(
                     tags$strong(
                       'Adjust heights of individual data tracks:'
                     )
                   ),
                   column(
                     6, 
                     wellPanel(
                       fluidRow(
                         column(
                           6,
                           sliderInput(
                             'expanded_plot_height_1', 'Methylation',
                             min = 1, max = 10, value = c(3), step = 1),
                         ),
                         column(
                           6,
                           sliderInput(
                             'expanded_plot_height_2', 'Transcripts',
                             min = 1, max = 10, value = c(1), step = 1)
                         )
                       )
                     )
                   )
                 ),
                 
                 
                 br(),
                 plotOutput('expanded', height = 900)),
        
        tabPanel("CpG Boxplots",
                 p('Show boxplots for specific CpGs:'),
                 
                 #cpgs search list related to the input
                 selectizeInput('cpg', 'CpGs related to the input',
                                choices = NULL,
                                multiple = TRUE),
                 
                 #action button to submit and reset input
                 actionButton('submit_cpg_gene_input', 'Submit'),
                 actionButton('reset_cpg_gene_input','Reset'),
                 
                 plotOutput('boxplot_selected_cpg')
                 
        )))),
  
  tags$h3('Annotation Table'),
  fluidRow(
    
    
    
    div(
      dataTableOutput('dt_anno_cpg', width = "98%"), 
      style = "font-size:85%",
      align = 'center')),
  
  #footer
  hr(),
  
  HTML(paste("[1] Only associated CpGs will be displayed when selecting a gene.This CpG-gene mapping is based on the Illumina-provided EPIC array annotation.",
             "[2] When selecting a gene, all CpGs that are mapped to any component of an associated transcript will be shown. Transcript components include 1-5Kb upstream of the TSS, the promoter (< 1Kb upstream of the TSS), 5’UTR, exons, introns, and 3’UTRs. Everything else are intergenic regions", 
             "[3] All DNA methylation data here is 850k array data. See primary publication for details on processing.",
             "[4] Annotations are in hg19 coordinates.",
             "[5] Differential methylation is defined as a bonferroni-adjusted p-value < 0.01, and an absolute difference in mean DNAm > 25%.",
             "[6] PMD regions coordinates are taken from D. I. Schroeder et al. 2013.",
             "[7] Differential methylation by cell type was conducted separated in first and third trimester samples for autosomal CpGs. Differentially methylated cytosines (DMCs) are defined as a bonferroni adjusted p-value < 0.01 and a mean difference in methylation > 10%. Differential methylation was conducted with a one-versus-all design.",
             
             
             sep = " <br> "))
)


## --------------------------------------------------------------------------
server <- function(input, output, session) {
  ############  searching lists ############
  #searching list of cpgs
  updateSelectizeInput(session = session, inputId = 'cpg_name', 
                       choices = cpg$value, server = TRUE)
  #searching list of genes
  updateSelectizeInput(session = session, inputId = 'gene_name', 
                       choices = gene$gene, selected = '',server = TRUE)
  #searching list of chromosomes
  chr <- unique(anno$chr) %>% str_sort(numeric = TRUE) %>% as.tibble()
  updateSelectizeInput(session = session, inputId = 'chr',
                       choices = chr$value, selected = '', server = TRUE)
  
  #get cpg sites related to the selected gene
  cpg_related_to_gene <- reactive(find_cpg_ids(get_gene_input()))
  cpg_related_to_region <- reactive(find_cpg_from_region(Chr = get_region_chr(), Start = get_region_start(), End = get_region_end()))
  #searching list (cpg sites related to selected gene)
  observe(
    { if(input$tabs == 'Gene Input')
    {updateSelectizeInput(session = session, 
                          inputId = 'cpg', 
                          choices = cpg_related_to_gene(), 
                          server = TRUE, 
                          selected = '')}
      else if(input$tabs == 'Region Input')
      {updateSelectizeInput(session = session, 
                            inputId = 'cpg', 
                            choices = cpg_related_to_region(), 
                            server = TRUE, 
                            selected = '')}
      else if(input$tabs == 'CpG Input')
      {updateSelectizeInput(session = session, 
                            inputId = 'cpg', 
                            choices = get_cpg_input(), 
                            server = TRUE, 
                            selected = '')}
      
      
    })
  
  ############  reset button ############  
  #reset the input and outputs when click the reset button
  observeEvent(input$reset_cpg_input, {reset('cpg_name')})
  observeEvent(input$reset_gene_input, {
    reset('gene_name') & reset('cpg_number_input') & reset('range_cpg') & reset('cpg')})
  observeEvent(input$reset_region_input, {
    reset('chr') & reset('start') & reset('end')& reset('cpg_number_input_region')})
  observeEvent(input$reset_gene_input, {input_value$input_value <- FALSE})
  observeEvent(input$reset_cpg_input, {input_value$input_value <- FALSE})
  observeEvent(input$reset_region_input, {input_value$input_value <- FALSE})
  #reset the input and output when switch between tabs
  observeEvent(input$tabs, {input_value$input_value <- FALSE})
  observeEvent(input$tabs, 
               reset('cpg_name') & reset('gene_name') & reset('chr') & reset('start') & 
                 reset('end') & reset('cpg_number_input') & reset('range_cpg') & reset('cpg'))
  
  ############  submit button ############
  #store inputs into reactive values
  input_value <- reactiveValues(input_value = FALSE)
  #get input values when click the submit button
  observeEvent(input$submit_cpg_input, {input_value$input_value <- input$submit_cpg_input})
  observeEvent(input$submit_gene_input, {input_value$input_value <- input$submit_gene_input})
  observeEvent(input$submit_gene_input, {
    reset('cpg_number_input')& reset('range_cpg') &reset('cpg')})
  observeEvent(input$submit_region_input, {input_value$input_value <- input$submit_region_input})
  observeEvent(input$submit_region_input, {reset('cpg_number_input_region') & reset('cpg')})
  
  
  ############  datatable output ############
  #make an empty datatable as default
  empty_df <- anno[0,] %>% rename_columns() %>% select(CpG:`Relation to Transcript Features`)
  #datatable to show cpg annotation information 
  
  output$dt_anno_cpg <- DT::renderDataTable({
    #output an empty datatable when no inputs
    if (input_value$input_value == FALSE)
      return(DT:: datatable(empty_df, 
                            options = list(paging = FALSE, 
                                           searching = FALSE)))
    #output a datatable when input a gene name
    isolate(if(input$tabs == 'Gene Input') {
      
      DT::datatable(
        make_datatable(gene_name = input$gene_name, 
                       cpg_name = NULL,
                       Chr = NULL, Start = NULL, End = NULL)[,input$checkbox,drop = FALSE],
        options = list(pageLength = 5),
        selection = 'single',
        rownames = FALSE) %>%
        formatStyle(1, cursor = 'pointer', color = 'blue')}
      
      #output a datatable when input cpgs
      else if(input$tabs == 'CpG Input') {
        
        DT::datatable(
          make_datatable(gene_name = NULL, cpg_name = input$cpg_name, 
                         Chr = NULL, Start = NULL, End = NULL)[,input$checkbox, drop = FALSE], 
          options = list(pageLength = 5),
          selection = 'single',
          rownames = FALSE) %>%
          formatStyle(1, cursor = 'pointer', color = 'blue')
      } else if(input$tabs == 'Region Input') {#output a datatable when input a region
        
        DT::datatable(
          make_datatable(
            gene_name = NULL, cpg_name = NULL,
            Chr = input$chr, Start = input$start, End = input$end
          )[,input$checkbox, drop = FALSE],
          options = list(pageLength = 5),
          selection = 'single',
          rownames = FALSE
        ) %>%
          formatStyle(1, cursor = 'pointer', color = 'blue')
        
      })
  }
  )
  
  ############  annoplot output: default to show first 250 cpgs############
  #gene cpg plots
  output$annoplot <- renderPlot(
    {
      #return empty when no inputs
      if (input_value$input_value  == FALSE) return()
      #output a plot when input a gene
      else if(input$tabs == 'Gene Input')
      {
        annoPlot_with_tracks(cpg_name = NULL, gene_name = input$gene_name, 
                             Chr = NULL, Start = NULL, End = NULL, 
                             Tissue_type = input$checkbox_tissue,
                             Sex_type = input$checkbox_sex,
                             Trimester_type = input$checkbox_trimester,
                             plot_heights = c(input$condensed_plot_height_1,
                                              input$condensed_plot_height_2,
                                              input$condensed_plot_height_3,
                                              input$condensed_plot_height_4))
      }
      #output a plot when input cpgs
      else if(input$tabs == 'CpG Input')
      {
        annoPlot_with_tracks(cpg_name = input$cpg_name, gene_name = NULL, 
                             Chr = NULL, Start = NULL, End = NULL, 
                             Tissue_type = input$checkbox_tissue,
                             Sex_type = input$checkbox_sex,
                             Trimester_type = input$checkbox_trimester,
                             plot_heights = c(input$condensed_plot_height_1,
                                              input$condensed_plot_height_2,
                                              input$condensed_plot_height_3,
                                              input$condensed_plot_height_4))
      }
      #output a plot when input a region
      else if(input$tabs == 'Region Input')
      {
        annoPlot_with_tracks(cpg_name =  NULL, gene_name = NULL, 
                             Start = input$start, End = input$end, Chr = input$chr, 
                             Tissue_type = input$checkbox_tissue,
                             Sex_type = input$checkbox_sex,
                             Trimester_type = input$checkbox_trimester,
                             plot_heights = c(input$condensed_plot_height_1,
                                              input$condensed_plot_height_2,
                                              input$condensed_plot_height_3,
                                              input$condensed_plot_height_4))
      }
      
      
    }
  )
  ############  annoplot output: choose cpg number (given a gene input) ############
  #get the input when click submit button
  get_gene_input <- eventReactive(input$submit_gene_input, {input$gene_name})
  get_cpg_input <- eventReactive(input$submit_cpg_input, {input$cpg_name})
  get_region_chr<- eventReactive(input$submit_region_input, {input$chr})
  get_region_start<- eventReactive(input$submit_region_input, {input$start})
  get_region_end<- eventReactive(input$submit_region_input, {input$end})
  
  #store cpg number into reactive value
  cpg_number <- reactiveValues(cpg_number = 250)
  #get the cpg number input
  observeEvent(input$submit_cpg_number, {cpg_number$cpg_number <- input$cpg_number_input})
  #reset cpg number into default(50) when select another gene
  observeEvent(input$submit_gene_input, {cpg_number$cpg_number = 250})
  #update the plot when change the cpg number to show
  observeEvent(input$submit_cpg_number, {output$annoplot <- renderPlot(
    if (input_value$input_value  == FALSE) return()
    else if(input$tabs == 'Gene Input')
    {
      annoPlot_with_tracks(cpg_name = NULL, gene_name = get_gene_input(), 
                           Chr = NULL, Start = NULL, End = NULL, 
                           cpg_number = cpg_number$cpg_number, 
                           Tissue_type = input$checkbox_tissue,
                           Sex_type = input$checkbox_sex,
                           Trimester_type = input$checkbox_trimester,
                           plot_heights = c(input$condensed_plot_height_1,
                                            input$condensed_plot_height_2,
                                            input$condensed_plot_height_3,
                                            input$condensed_plot_height_4))
      
    }
    #output a plot when input cpgs
    else if(input$tabs == 'CpG Input')
    {
      annoPlot_with_tracks(cpg_name = get_cpg_input(), gene_name = NULL, 
                           Chr = NULL, Start = NULL, End = NULL, 
                           Tissue_type = input$checkbox_tissue,
                           Sex_type = input$checkbox_sex,
                           Trimester_type = input$checkbox_trimester,
                           plot_heights = c(input$condensed_plot_height_1,
                                            input$condensed_plot_height_2,
                                            input$condensed_plot_height_3,
                                            input$condensed_plot_height_4))
    }
    #output a plot when input a region
    else if(input$tabs == 'Region Input')
    {
      annoPlot_with_tracks(cpg_name =  NULL, gene_name = NULL, 
                           Start = get_region_start(), 
                           End = get_region_end(), 
                           Chr = get_region_chr(), 
                           Tissue_type = input$checkbox_tissue,
                           Sex_type = input$checkbox_sex,
                           Trimester_type = input$checkbox_trimester,
                           plot_heights = c(input$condensed_plot_height_1,
                                            input$condensed_plot_height_2,
                                            input$condensed_plot_height_3,
                                            input$condensed_plot_height_4))
    }
  )
  })
  ############  annoplot output: choose cpg number (given a region input) ############
  #store cpg number into reactive value
  cpg_number_region <- reactiveValues(cpg_number = 250)
  #get the cpg number input
  observeEvent(input$submit_cpg_number_region, {cpg_number_region$cpg_number <-
    input$cpg_number_input_region})
  #reset cpg number into default(50) when select another gene
  observeEvent(input$submit_region_input, {cpg_number_region$cpg_number = 250})
  #update the plot when change the cpg number to show
  observeEvent(input$submit_cpg_number_region, {output$annoplot <- renderPlot(
    if (input_value$input_value  == FALSE) return()
    else if(input$tabs == 'Gene Input')
    {
      annoPlot_with_tracks(cpg_name = NULL, gene_name = get_gene_input(),
                           Chr = NULL, Start = NULL, End = NULL, 
                           Tissue_type = input$checkbox_tissue,
                           Sex_type = input$checkbox_sex,
                           Trimester_type = input$checkbox_trimester,
                           plot_heights = c(input$condensed_plot_height_1,
                                            input$condensed_plot_height_2,
                                            input$condensed_plot_height_3,
                                            input$condensed_plot_height_4))
      
    }
    #output a plot when input cpgs
    else if(input$tabs == 'CpG Input')
    {
      annoPlot_with_tracks(cpg_name = get_cpg_input(), gene_name = NULL, 
                           Chr = NULL, Start = NULL, End = NULL, 
                           Tissue_type = input$checkbox_tissue,
                           Sex_type = input$checkbox_sex,
                           Trimester_type = input$checkbox_trimester,
                           plot_heights = c(input$condensed_plot_height_1,
                                            input$condensed_plot_height_2,
                                            input$condensed_plot_height_3,
                                            input$condensed_plot_height_4))
    }
    #output a plot when input a region
    else if(input$tabs == 'Region Input')
    {
      annoPlot_with_tracks(cpg_name =  NULL, gene_name = NULL, 
                           Start = get_region_start(), 
                           End = get_region_end(), 
                           Chr = get_region_chr(), 
                           cpg_number = cpg_number_region$cpg_number, 
                           Tissue_type = input$checkbox_tissue,
                           Sex_type = input$checkbox_sex,
                           Trimester_type = input$checkbox_trimester,
                           plot_heights = c(input$condensed_plot_height_1,
                                            input$condensed_plot_height_2,
                                            input$condensed_plot_height_3,
                                            input$condensed_plot_height_4))
    }
  )
  })
  
  ############  annoplot output: choose a set of cpg (given a gene input)############
  #store cpg_start and cpg_end into reactive value
  cpg_start <- reactiveValues(cpg_start = 1)
  cpg_end <- reactiveValues(cpg_start = 500)
  #get the cpg set input
  observeEvent(input$submit_cpg_set, {cpg_start$cpg_start <- input$range_cpg[1]
  cpg_end$cpg_end <- input$range_cpg[2]})
  #reset cpg set when select another gene
  observeEvent(input$submit_gene_input, {cpg_start$cpg_start = 1
  cpg_end$cpg_end = 500})
  #update the plot when change the cpg number to show
  observeEvent(input$submit_cpg_set, {output$annoplot <- renderPlot(
    if (input_value$input_value  == FALSE) return()
    else if(input$tabs == 'Gene Input')
    {
      annoPlot_with_tracks(cpg_name = NULL, gene_name = get_gene_input(), 
                           Chr = NULL, Start = NULL, End = NULL, 
                           first_cpg = cpg_start$cpg_start, end_cpg = cpg_end$cpg_end, 
                           Tissue_type = input$checkbox_tissue,
                           Sex_type = input$checkbox_sex,
                           Trimester_type = input$checkbox_trimester,
                           plot_heights = c(input$condensed_plot_height_1,
                                            input$condensed_plot_height_2,
                                            input$condensed_plot_height_3,
                                            input$condensed_plot_height_4))
      
    }
    #output a plot when input cpgs
    else if(input$tabs == 'CpG Input')
    {
      annoPlot_with_tracks(cpg_name = get_cpg_input(), gene_name = NULL, 
                           Chr = NULL, Start = NULL, End = NULL, 
                           Tissue_type = input$checkbox_tissue,
                           Sex_type = input$checkbox_sex,
                           Trimester_type = input$checkbox_trimester,
                           plot_heights = c(input$condensed_plot_height_1,
                                            input$condensed_plot_height_2,
                                            input$condensed_plot_height_3,
                                            input$condensed_plot_height_4))
    }
    #output a plot when input a region
    else if(input$tabs == 'Region Input')
    {
      annoPlot_with_tracks(cpg_name =  NULL, gene_name = NULL, 
                           Start = get_region_start(), 
                           End = get_region_end(), 
                           Chr = get_region_chr(), 
                           first_cpg = cpg_start$cpg_start, 
                           end_cpg = cpg_end$cpg_end, 
                           Tissue_type = input$checkbox_tissue,
                           Sex_type = input$checkbox_sex,
                           Trimester_type = input$checkbox_trimester,
                           plot_heights = c(input$condensed_plot_height_1,
                                            input$condensed_plot_height_2,
                                            input$condensed_plot_height_3,
                                            input$condensed_plot_height_4))
    }
  )
  })
  
  
  
  ############  Expanded view output ############
  ############  (given a gene input) ############
  #gene cpg plots 
  output$expanded <- renderPlot(
    {
      #return empty when no inputs
      if (input_value$input_value  == FALSE) return()
      
      #output a plot when input a gene
      else if(input$tabs == 'Gene Input')
      {
        expanded_view(gene_name = input$gene_name, 
                      
                      start_extend = 2500, # how far to extend the plot left, in genomic coordinates 
                      end_extend = 2500,   # how far to extend right
                      show_transcript_gene = TRUE, # show transcripts associated gene symbol
                      breaks_width = 15000, # determines the spacing of x axis breaks
                      trimester = input$checkbox_trimester, # can be Third, First, or Both
                      sex = input$checkbox_sex, # can be Male, Female, or Both
                      tissue = input$checkbox_tissue,
                      plot_sizes = c(input$expanded_plot_height_1,
                                     input$expanded_plot_height_2))
      }
    }
  )
  
  ############  boxplot output ############
  ##retrive cpg name when click the datatable
  get_cpg_from_dt <- eventReactive(input$submit_cpg_gene_input, {input$cpg})
  
  #produce a boxplot based on the selected cpg from datatable
  observeEvent(input$dt_anno_cpg_cell_clicked, 
               {
                 c_info = input$dt_anno_cpg_cell_clicked
                 
                 output$boxplot_selected_cpg <- renderPlot(
                   {
                     if (!is.null(c_info$value))
                     {
                       boxplot_selected(c_info$value, sex_type = input$checkbox_sex)
                     }
                     else if(!is.null(get_cpg_from_dt()))
                     {
                       boxplot_selected(get_cpg_from_dt(), sex_type = input$checkbox_sex)
                     }
                   }
                 )
               })
  #produce a boxplot based on the selected cpg using search box
  observeEvent(input$submit_cpg_gene_input, {output$boxplot_selected_cpg <- renderPlot(
    {
      boxplot_selected(get_cpg_from_dt(), sex_type = input$checkbox_sex)
    }
  )})
  
  ############  reset boxplot output ############
  #reset the associated cpg inputs in search box
  observeEvent(input$reset_cpg_gene_input, {reset('cpg')})
  #reset the boxplot when click the reset button
  observeEvent(input$reset_cpg_gene_input, {output$boxplot_selected_cpg <- renderPlot(
    {
      return()
    }
  )})
  #reset the boxplot when switch between tabs
  observeEvent(input$tabs, {output$boxplot_selected_cpg <- renderPlot(
    {
      return()
    }
  )})
  #reset the boxplot when change into another gene
  observeEvent(input$submit_gene_input, {output$boxplot_selected_cpg <- renderPlot(
    {
      return()
    }
  )})
  #reset the boxplot when reset gene input
  observeEvent(input$reset_gene_input, {output$boxplot_selected_cpg <- renderPlot(
    {
      return()
    }
  )})
  
  
}


## --------------------------------------------------------------------------
shinyApp(ui, server)

