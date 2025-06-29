rename_columns <- function(x) {
  dplyr::rename(
    x,
    
    # columns to rename
    CpG = cpg,
    Chr = chr,
    Pos = pos,
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

#produce a datatable which shows cpg annotation
make_datatable <- function(cpg_name = NULL, gene_name = NULL, 
                           Chr = NULL, Start = NULL, End = NULL){
  if (!is.null(cpg_name)) { #get cpg annotation 
    anno_cpg <- filter(anno, cpg %in% cpg_name) %>% arrange(desc(pos))
  }
  else if (!is.null(gene_name)) {
    #get cpg sites related to input genes
    cpg_name <- find_cpg_ids(gene_name)
    
    #get cpg annotation 
    anno_cpg <- filter(anno, cpg %in% cpg_name)%>% arrange(desc(pos))
  }
  
  else if (!is.null(Chr) & !is.null(Start) & !is.null(End))
  {
    #get cpg sites within a given region
    cpg_name <- find_cpg_from_region(Chr, Start, End)
    
    #get cpg annotation 
    anno_cpg <- filter(anno, cpg %in% cpg_name)%>% arrange(desc(pos)) 
  }
  
  anno_cpg <- anno_cpg %>%
    rename_columns()
}

annoPlot_with_tracks <- function(cpg_name = NULL, 
                                 gene_name = NULL, 
                                 Chr = NULL, 
                                 Start = NULL, 
                                 End = NULL, 
                                 cpg_number = 500, 
                                 first_cpg = NULL, 
                                 end_cpg = NULL, 
                                 Tissue_type = unique(pDat$Tissue),
                                 Sex_type = c('F', 'M'),
                                 Trimester_type = 'Third',
                                 plot_heights = c(10, 5, 5, 5)){
  
  ######## FILTER ANNOTATION TO RELEVANT CPGS, DEPENDING ON INPUT #################
  if (!is.null(gene_name) & is.null(first_cpg) & is.null(end_cpg)){
    
    # If a gene symbol is supplied,
    # filter to gene in annotation
    anno_gene <- anno %>%
      filter(grepl(paste0('(?<!-)\\b', gene_name, '\\b(?!-)'), genes_symbol, perl = TRUE),
             cpg %in% rownames(betas)) %>%
      arrange(desc(pos)) %>% slice(1:cpg_number)
    
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
      arrange(desc(pos)) %>% slice(first_cpg: end_cpg)
    
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
      arrange(chr, desc(pos))
    
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
      arrange(chr, desc(pos)) %>% slice(1:cpg_number)
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
      arrange(chr, desc(pos)) %>% slice(first_cpg : end_cpg)
    
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
    scales::comma(min(anno_gene$pos)), '-',
    scales::comma(max(anno_gene$pos))
  )
  
  if (is.null(cpg_name)){
    title <- paste0(gene_name, ', ', title)
  } 
  
  
  #
  # now we wrangle the betas/annotations/pdata information into a format usable for plotting
  #
  betas_gene <- t(betas[anno_gene$cpg,,drop = FALSE]*100) %>% 
    as_tibble() %>% 
    mutate(id = colnames(betas)) %>%
    
    # reshape into longer format
    pivot_longer(cols = -id,
                 names_to = 'cpg', 
                 values_to = 'beta')  %>%
    
    # separate out mean and sd
    separate(id, into = c('var', 'Trimester', 'Tissue'),
             sep = '_') %>%
    
    # put sd and mean into columns
    pivot_wider(
      names_from = var,
      values_from = beta) %>%
    
    # apply filtering criteria
    filter(Trimester %in% Trimester_type) %>%
    
    
    #calculate global mean for each cpg
    group_by(Trimester, cpg) %>%
    mutate(mean_cpg = mean(mean)) %>%
    
    # calculate lower, upper, and mean difference
    ungroup() %>%
    mutate(lower = mean-sd, upper = mean+sd,
           mean_diff = mean-mean_cpg) %>%
    
    # add cpg info
    left_join(anno_gene, by = 'cpg') %>%
    
    # order on cpg position
    ungroup() %>%
    arrange(desc(pos)) %>%
    mutate(cpg = factor(cpg, levels = unique(cpg)),
           Trimester_Tissue = paste0(Trimester, ' - ', Tissue)) %>%
    filter(Tissue %in% Tissue_type)
  
  ##### GENERATE INDIVIDUAL PLOT TRACKS #####
  # common theme settings across each plot
  
  font_size <- 14
  point_size <- 3
  line_size <- 1.25
  
  theme_custom <- function(base_size = font_size){
    theme_bw(base_size = base_size) %+replace%
      theme(
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = font_size-2, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0, size = 16),
        
        legend.text = element_text(size = font_size),
        legend.box = "horizontal",
        legend.position = 'right',
        legend.justification = 'left',
        
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        
        panel.spacing.y = unit(0.1, 'cm'),
        
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold', hjust = 0,
                                  margin = margin(b = 4)),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  }
  
  ########## Generate  methylation track
  # Term
  p1a <- betas_gene %>%
    filter(Trimester == Trimester_type) %>%
    
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
    
    # add trimester subtitle
    facet_wrap(vars(Trimester),
               ncol = 1, 
               labeller = labeller(Trimester = c('First' = 'DNA methylation, First Trimester', 
                                                 'Third' = 'DNA methylation, Term')),
               scales = 'free') +
    
    # cosmetics
    theme_custom() +
    scale_y_continuous(limits = c(0, 100), 
                       expand = c(0, 0), 
                       breaks = c(0, 50, 100),
                       labels = function(x)paste0(x, '%')) +
    scale_x_discrete(expand = c(0,0.5))+
    scale_color_manual(values= color_code_tissue[unique(betas_gene$Tissue)],
                       guide = guide_legend(ncol = 1, override.aes = list(size = 4)),
                       labels = function(x)gsub(' cs', '', x)) +
    labs(y = '', 
         x = '',
         title = title,
         color = '')
  
  # difference fro mean
  p1b <- betas_gene %>%
    filter(Trimester == Trimester_type) %>%
    ggplot(aes(x = cpg, color = Tissue, fill = Tissue, y = mean_diff)) +
    
    geom_hline(yintercept = 0, color = 'black', size = line_size*0.75) +
    geom_point(size = point_size/2, show.legend = FALSE) +
    geom_line(aes(group = Tissue)) +
    facet_wrap(vars(Trimester),
               ncol = 1, 
               labeller = labeller(Trimester = c('First' = 'Difference from CpG mean', 
                                                 'Third' = 'Difference from CpG mean'))) +
    theme_custom() +
    theme(panel.grid.minor.y = element_blank()) +
    scale_color_manual(values= color_code_tissue[unique(betas_gene$Tissue)],
                       guide = guide_none(),
                       labels = function(x)gsub(' cs', '', x)) +
    scale_y_continuous(labels = function(x)paste0(x, '%'),
                       breaks = c(-100, -50, 0, 50, 100),
                       limits = c(-100, 100))+
    labs(x = '', y = '') 
  
  # plot cpg annotations
  p2 <- betas_gene %>% 
    
    # First we need to process the annotation data
    mutate(# absent/presence for different genomic elements
      enhancer = !is.na(enhancers_id),
      pmd = !is.na(pmd_id)) %>%
    mutate_if(is.logical, as.character) %>%
    select(cpg, enhancer, pmd, 
           imprinted_dmr_general, imprinted_dmr_placenta, cpg_id, pos) %>%
    
    # reshape
    pivot_longer(cols = c(-cpg, -pos),
                 names_to = 'cpg_element', 
                 values_to = 'presence') %>%
    distinct() %>%
    
    # arrange plot orde for each track
    arrange(desc(pos)) %>%
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
    theme_custom() +
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
    select(cpg, pos, genes_id, genes_tx_id, genes_symbol) %>%
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
    mutate(genes_symbol_end = min(pos)) %>%
    ungroup() %>%
    
    # reorder gene symbols by their earliest end position
    mutate(genes_symbol = fct_reorder(genes_symbol, genes_symbol_end)) %>%
    
    # reorder transcripts based on gene symbol levels and then by their end position
    arrange(genes_symbol, pos) %>%
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
    theme_custom()+
    theme(axis.text.x = element_text(angle = 90)) +
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
  if (nrow(anno_gene) > 100) {
    p3 <- p3 + theme(axis.text.x = element_blank())
    
  }
  
  # combine tracks plot
  (p1a + plot_layout(guides = 'keep')) / 
    p1b / 
    p2 / 
    p3 + plot_layout(heights =c(plot_heights))
}

# data wrangling -------
find_cpg_ids <- function(gene_name){
  cpg_name <- anno %>%
    filter(cpg %in% rownames(betas),
           grepl(paste0('(?<!-)\\b', gene_name, '\\b(?!-)'), genes_symbol, perl = TRUE)) %>%
    pull(cpg) %>% 
    unique()
  
  cpg_name
}

find_cpg_from_region <- function(Chr, Start, End){
  
  cpg_name <- anno %>% 
    filter(chr == Chr, pos > Start, pos < End) %>%
    pull(cpg) %>%
    unique()
  
  cpg_name <- cpg_name[cpg_name %in% rownames(betas)]
  
  cpg_name
}

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
    arrange(desc(pos)) %>% 
    mutate(cpg = factor(cpg, levels = unique(cpg)))
  cpg
}

# vis ------
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

expanded_view <- function(gene_name = NULL, 
                          
                          start_extend = 2500, # how far to extend the plot left, in genomic coordinates 
                          end_extend = 2500,   # how far to extend right
                          
                          show_transcript_gene = TRUE, # show transcripts associated gene symbol
                          breaks_width = 15000, # determines the spacing of x axis breaks
                          
                          highlight_region = NULL, # a df with start end columns, indicating regions to highlight
                          
                          trimester = 'Third', # can be Third, First, or Both
                          sex = c('F', 'M'), # can be Male, Female, or Both
                          tissue = c("Endothelial cs", "Villi", "Hofbauer cs", 
                                     "Trophoblasts cs", "Stromal cs"),
                          plot_sizes = c(3, 1) # a 2-element numeric vector indicating the relative plot heights 
                          # of the methylation plot to the transcript plot
                          
                          # required components:
                          # pDat with columns Sample_Name, Trimester, Sex, and Tissue
                          # betas_filt data.frame with column names == pDat$Sample_Name
                          # anno_annotatr data frame that is contains a column cpg, genes_symbol, chr, start, end
                          # the cpgs in anno_annotatr and betas_filt need to be exactly the same, and same order
                          # annotation data.frame that I built from annotatr, and slightly cleaned
                          
) {
  
  # subset to cpgs associated with gene
  anno_gene_subset <- anno %>%
    filter(grepl(paste0('\\<', gene_name, '\\>'), genes_symbol)) %>% 
    dplyr::select(chr, pos) 
  
  c <- unique(anno_gene_subset$chr)
  s <- min(anno_gene_subset$pos) # first cpg
  e <- max(anno_gene_subset$pos) # last cpg
  s_extend <- s - start_extend
  e_extend <- e + end_extend
  
  region_of_interest <- anno %>%
    filter(cpg %in% rownames(betas),
           chr == c, pos > s_extend, pos < e_extend)
  
  # subset betas
  b <- betas[region_of_interest$cpg,] %>%
    as.data.frame() %>%
    
    # make long
    bind_cols(cpg = rownames(.), .) %>%
    pivot_longer(cols = -cpg,
                 names_to = 'id',
                 values_to = 'beta') %>%
    
    # separate out mean and sd
    separate(id, into = c('var', 'Trimester', 'Tissue'),
             sep = '_') %>%
    
    # put sd and mean into columns
    pivot_wider(
                names_from = var,
                values_from = beta) %>%
    
    # apply filtering criteria
    filter(Trimester %in% trimester,
           Tissue %in% tissue) %>%
    
    
    left_join(anno, by = 'cpg')
  
  font_size <- 14
  point_size <- 3
  line_size <- 1.25
  
  # methylation plot
  p1 <- b %>%
    
    ggplot(aes(x = pos)) +
    geom_linerange(ymin = 0, aes(ymax = mean, color = Tissue)) +
    
    facet_wrap(vars(Tissue), ncol = 1, strip.position = 'left',
               labeller = labeller(Tissue = function(x)gsub(' cs', '', x))) +
    
    theme_bw(base_size = font_size) +
    theme(strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1),
          strip.placement = 'outside',
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_line(size = 0.5),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0, 'lines'),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
          #panel.grid = element_blank()
    ) +
    scale_y_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0,0), labels = NULL) +
    scale_x_continuous(limits = c(s_extend, e_extend), expand = c(0,0), 
                       #breaks = scales::breaks_width(15000),
                       labels = scales::number) +
    scale_color_manual(values = colors, guide = 'none') +
    labs(y = '', x = '', title = gene_name, 
         subtitle = paste0(c, ':', scales::comma(s), '-', scales::comma(e)))
  
  p2 <- anno_annotatr %>%
    filter(seqnames == as.character(c), start > s_extend, end < e_extend,
           type %in% c('exons', 'introns', 'promoters')) %>%
    
    mutate(symbol = ifelse(is.na(symbol), '', symbol)) %>%
    ggplot(aes(x = start, xend = end)) +
    geom_segment(y= 0.5, yend = 0.5, aes(size = type, color = type))  +
    facet_wrap(~symbol + tx_id, ncol = 1, strip.position = 'left') +
    theme_bw(base_size = font_size)+
    theme(strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1),
          strip.placement = 'outside',
          legend.position = 'bottom',
          plot.title =  element_blank(),
          plot.subtitle = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(size = 0.5),
          panel.spacing = unit(0, 'lines'),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.border = element_blank(),
          legend.key.size = unit(0.5, 'cm'),
          legend.spacing = unit(0.25, 'cm'),
          legend.direction = 'vertical',
          legend.margin = margin(t = -20,0,0,0)) +
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
                       #breaks = scales::breaks_width(15000)
    ) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = 'position', size = '', color = '')
  
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