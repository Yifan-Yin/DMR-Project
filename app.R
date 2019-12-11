#' ---
#' title: "DMC_alltest Visualization"
#' output: html_document
#' editor_options: 
#'   chunk_output_type: console
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' #package setup
#' 
## ------------------------------------------------------------------------------
library(tidyverse)
library(egg)
library(here)

#' 
#' # Setup
#' 
#' ## Load data
#' 
## ------------------------------------------------------------------------------
# sample data
pDat <- readRDS(here::here('data', '3_1_pDat_filt.rds'))

# methylation data
betas <- readRDS(here::here('data', '3_1_betas_filt.rds'))

# annotation
anno <- readRDS(here::here('data', 'annotation.rds'))

#color code
color_code <- readRDS(here::here('data', '2_3_color_code.rds'))
color_code_tissue <- setNames(color_code$Colors_Tissue, color_code$label)
colors <- color_code_tissue[unique(pDat$Tissue)]

#' 
#' ## Rename column names
#' 
#' For making displayed column names more sensible.
#' 
#' This is a function that renames column names. I apply it to the output of several functions that generate the annotation table.
#' 
## ------------------------------------------------------------------------------
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

#' 
#' ## Column names
#' 
#' This is for checkboxGroupInput
#' 
## ------------------------------------------------------------------------------
#column names 
column_names <- anno %>% 
  select(cpg, chr, start, end, 
         cpg_id, enhancers_id, genes_symbol, genes_tx_id, genes_id, pmd_id, 
         imprinted_dmr_general, imprinted_dmr_placenta) %>%
  rename_columns() %>%
  names()

#' 
#' ## Generate gene/cpg search list
#' 
## ------------------------------------------------------------------------------
#list of genes for search list
gene <- tibble(gene = str_split(anno$genes_symbol, ', ') %>%
                 unlist() %>%
                 unique() %>%
                 na.omit() %>%
                 sort())

#pull out cpg name for search list use
cpg <- unique(anno$cpg) 
cpg <- cpg[cpg %in% rownames(betas)] %>% as_tibble()

#' 
#' # 1. datatable which shows cpg annotation
#' 
#' ## Main function
#' 
## ------------------------------------------------------------------------------
#datatable which shows cpg annotation (when inputs have both cpgs and genes)
cpg_info_total <- function(cpg_name = NULL, gene_name = NULL){
  
  
  if (is.null(cpg_name)){
    
    # when input is gene symbol
    cpgs <- fn_cpg(gene_name)
    info <- filter(anno, cpg %in% cpgs)%>% arrange(desc(start))
  }
  else if (gene_name == ''){
    
    # when input are cpg ids
    info <- filter(anno, cpg %in% cpg_name) %>% arrange(desc(start))
  }
  else {
    stop('please enter either cpgs or a gene')
  }
  
  info %>%
    rename_columns()
}

#' 
#' # 2. Main gene/cpg plot with tracks
#' 
#' This is the middle plot with all the annotation tracks.
#' 
## ------------------------------------------------------------------------------
stat_plot <- function(cpg_name = NULL, gene_name = NULL){
  
  ######## FILTER ANNOTATION TO RELEVANT CPGS, DEPENDING ON INPUT #################
  if (is.null(cpg_name)){
    
    # If a gene symbol is supplied,
    # filter to gene in annotation
    anno_gene <- anno %>%
      filter(grepl(paste0('(?<!-)\\b', gene_name, '\\b(?!-)'), genes_symbol, perl = TRUE),
             cpg %in% rownames(betas)) %>%
      arrange(desc(start))
    
    # stop if gene symbol is invalid
    if (nrow(anno_gene) == 0) {
      stop('No cpgs found, maybe they were filtered out.')
    }
  } else if (gene_name == ''){
    
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
    
  } else {
    stop('Enter valid cpg(s) OR one specific gene.')
  }
  
  ######## SETUP INPUT TO PLOTS ###############
  #
  # Title
  #
  title <- paste0(
    unique(anno_gene$chr), ':',
    min(anno_gene$start), '-',
    max(anno_gene$start)
  )
  
  if (is.null(cpg_name)){
    title <- paste0(gene_name, ', ', title)
  } 
  
  
  #
  # now we wrangle the betas/annotations/pdata information into a format usable for plotting
  #
  betas_gene <- t(betas[anno_gene$cpg,,drop = FALSE]*100) %>% as_tibble %>% 
    mutate(Sample_Name = colnames(betas)) %>%
    
    # reshape into longer format
    pivot_longer(cols = -Sample_Name,
                 names_to = 'cpg', 
                 values_to = 'beta')  %>%
    
    # add tissue and trimester info
    left_join(pDat %>% select(Sample_Name, Trimester, Tissue), by = 'Sample_Name') %>%
    
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
           Trimester_Tissue = paste0(Trimester, ' - ', Tissue))
  
  ##### GENERATE INDIVIDUAL PLOT TRACKS #####
  font_size <- 12
  point_size <- 3
  line_size <- 1.25
  
  # Generate  methylation track
  p1 <- betas_gene %>%
    ggplot() +
    
    # geom_ribbon is to generate the alternating shaded background
    geom_ribbon(aes(x = as.numeric(cpg), ymin = 0, ymax = 25),
                fill = 'grey', alpha = 0.15)+
    geom_ribbon(aes(x = as.numeric(cpg), ymin = 50, ymax = 75),
                fill = 'grey', alpha = 0.15)+
    
    # geom_linerange and geom_point is for the actual methylation data
    geom_linerange(alpha = 0.5, size = line_size, 
                   aes(x = cpg, ymin =lower, ymax = upper, color = Tissue),
                   show.legend = FALSE) +
    geom_point(alpha = 0.75, aes(x = cpg, y = mean, color = Tissue), size = point_size) +
    
    # we want to facet by trimester
    facet_wrap(vars(Trimester), ncol = 1, labeller = label_both) +
    
    # cosmetics
    theme_bw(base_size = font_size) +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0, size = 16),
          
          legend.text = element_text(size = 14),
          
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          
          panel.spacing.y = unit(0.5, 'cm'),
          
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold', hjust = 0),
          
          plot.margin=margin(l=-0.1, unit="cm")) +
    scale_y_continuous(limits = c(0, 100), 
                       expand = c(0, 0), 
                       labels = function(x)paste0(x, '%')) +
    scale_x_discrete(expand = c(0.01,0.01))+
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
    mutate(facet_title = 'Differential methylation') %>%
    separate(Group, c('Trimester', 'Tissue'), '_') %>%
    
    # convert axes to numeric for geom_path to work
    mutate(cpg_num = as.numeric(cpg),
           Tissue_num = as.numeric(as.factor(Tissue)),
           
           # Add tissue-specific color when significant
           Significant_num = if_else(Significant, Tissue_num, NA_real_)) %>% 
    
    {
      ggplot(data = ., aes(x = cpg_num, y = Significant_num, color = Tissue)) +
        geom_point(data = . %>% 
                     filter(!is.na(Significant_num)), size = point_size) +
        geom_path(na.rm = TRUE, size = line_size-0.25, alpha = 1) +
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
              
              plot.margin = margin(l=-0.1,unit="cm"),
              plot.background = element_rect(fill = '#f7f7f7'),
              
              strip.background = element_blank(),
              strip.text.y = element_text(face = 'bold', vjust = 0.5, hjust = 0.5, angle = 180),
              strip.text.x = element_text(face = 'bold', hjust = 0),
              strip.placement = 'outside') +
        scale_color_manual(values = color_code_tissue[unique(betas_gene$Tissue)],
                           na.value = '#f7f7f7',
                           guide = 'none') +
        scale_y_continuous(expand = c(0.15,0.15),
                           breaks = c(1, 2, 3, 4, 5), 
                           labels = betas_gene$Tissue %>% 
                             as.factor() %>% 
                             levels() %>%
                             gsub(' cs', '', .)) +
        labs(x = '', y = '')
      
    }
  
  # plot cpg annotations
  p2 <- betas_gene %>% 
    
    # First we need to process the annotation data
    mutate(# absent/presence for different genomic elements
      enhancer = !is.na(enhancers_id),
      pmd = !is.na(pmd_id),
      imprinted_dmr_general = !is.na(imprinted_dmr_general),
      imprinted_dmr_placenta = !is.na(imprinted_dmr_placenta)) %>%
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
          
          plot.margin = margin(l=-0.1,unit="cm"),
          
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
    scale_x_discrete(expand = c(0,0)) +
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
          
          plot.margin = margin(l=-0.1,unit="cm"),
          
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
    scale_x_discrete(expand = c(0,0)) +
    coord_cartesian(clip = "off") +
    labs(x = '', y = '', fill = '')
  
  # the size of each plot depends on the number of unique transcripts
  # More transcripts means a longer transcript track
  # so the plot height depends on the number of transcripts
  
  # number of transcripts
  n_trans <- betas_gene$genes_tx_id %>% str_split(', ') %>% unlist() %>% unique %>% length()
  
  # size of p1
  p1_h <- 40
  
  # combine tracks plot
  ggarrange(p1, p_dmc, p2, p3, ncol = 1, heights = c(p1_h-(n_trans-6), 9, 6, n_trans))
}


#' 
#' # 3. CpG boxplots
#' 
#' These functions define the individual CpG plots when the user selects specific CpGs.
#' 
#' ## Data wrangling functions
#' 
#' These set of functions help extract the relevant information given user input, and then manipulates the data into a format usable by the plotting functions.
#' 
#' ### Given a gene symbol, return associated cpg IDs
#' 
## ------------------------------------------------------------------------------
fn_cpg <- function(gene_name){
  cpg_name <- anno %>%
    filter(cpg %in% rownames(betas),
           grepl(paste0('(?<!-)\\b', gene_name, '\\b(?!-)'), genes_symbol, perl = TRUE)) %>%
    pull(cpg) %>% 
    unique()
  
  cpg_name
}

#' 
#' ### Given cpgs, pull out DNA methylation from betas dataframe
#' 
#' Given a cpg, pull the associated DNA methylation at these cpgs, along with tissue and trimester information.
#' 
## ------------------------------------------------------------------------------
sample_info_cpg <- function(cpg_name){
  
  cpg <- betas[cpg_name,]
  
  #when have only one cpg input
  if (length(cpg_name) ==1){
    cpg <- as.data.frame(cpg) %>% t()
    cpg <- cpg[,colnames(cpg) %in% pDat$Sample_Name]
    cpg <- as.data.frame(cpg)
    colnames(cpg) <- cpg_name
    cpg <- cpg %>% cbind(pDat) %>% gather(cpg, betas, contains('cg'))
  }
  
  #when have more than one cpg inputs
  else{
    cpg <- cpg[,colnames(cpg) %in% pDat$Sample_Name]
    cpg <- cpg %>% t() %>% cbind(pDat) %>% gather(cpg, betas, contains('cg'))
  }
  cpg <- cpg %>% 
    left_join(anno, by = 'cpg') %>% 
    arrange(desc(start)) %>% 
    mutate(cpg = factor(cpg, levels = unique(cpg)))
  cpg
}

#' 
#' Given a gene, pull the relevant cpgs, then betas + tissue/trimester information
#' 
## ------------------------------------------------------------------------------
sample_info_gene <- function(gene_name){
  #get related cpg sites
  cpg_name <- fn_cpg(gene_name)
  #get sample cell types and tri
  cpg_info <- sample_info_cpg(cpg_name) 
  cpg_info
}

#' 
#' ## Plot function
#' 
## ------------------------------------------------------------------------------
#boxplots which show beta values based on cell types
boxplot_individual<- function(cpg_name){
  
  #given cpg, get sample betas, cell type  and trimester
  cpg <- sample_info_cpg(cpg_name)
  
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

#' ---
#' title: "shiny_DMR"
#' output: html_document
#' ---
#' 
#' #set up
## ------------------------------------------------------------------------------
#install.packages('rsconnect')
#install.packages('shinyjs')
library(rsconnect)
library(shiny)
rsconnect::setAccountInfo(name='yifan7', token='52E1037550A7AA68252579A98EE8AA38', secret='+NlrJysSVZHa5LQQb9slmYRDXExPfjF0QbxMjUH1')
library(DT)
library(shinyjs)

#' 
#' #shiny
## ------------------------------------------------------------------------------
ui <- fluidPage(
  theme = shinythemes::shinytheme("lumen"),
  
  useShinyjs(),
  #main title
  titlePanel("Placental Cell Methylome Browser"),
  hr(), #create a line
  
  #UI layout
  sidebarLayout(
    sidebarPanel(
      
      # relative width out of 12 units
      width = 3,
      
      #app description
      p('Welcome to the Placental Cell Methylome Browser!'),
      p('This app allows visualization of DNA methylation across various placental tissues and cell types.'),
      p('To start, select CpGs or a gene.'),
      hr(),
      
      #cpg name input
      selectizeInput('name', 'CpG Sites', 
                     choices = NULL,
                     multiple = TRUE,
                     options = list(placeholder='Enter CpG site names' )),
      #gene input  
      selectizeInput('gene', 'Gene', 
                     choices = NULL, 
                     multiple = FALSE,
                     options = list(placeholder='Enter a gene name',
                                    onInitialize = I('function() { this.setValue(""); }'))),
      
      actionButton("submit", "Submit"),
      actionButton('reset','Reset'),
      br(), # add a space line
      hr(), # create a horizontal line
      
      p('If a specific gene is selected, then you can also display specific associated-CpGs'),
      
      #cpgs related to selected gene
      selectizeInput('cpg', 'CpGs related to the selected gene', 
                     choices = NULL,
                     multiple = TRUE),
      actionButton('submit2', 'Submit'),
      actionButton('reset2','Reset'),
      
      br(),
      hr(),#create a horizontal line
      checkboxGroupInput('checkbox', 
                         'Columns to show in CpG annotation table', 
                         choices = column_names, 
                         selected = column_names[-10:-12])
    ), 
    
    mainPanel(
      hr(),
      h6("CpG Information. Click on a CpG ID to show a boxplot at the bottom."),
      
      #create a data table to show cpg site information 
      div(dataTableOutput('info'), style = "font-size:85%"),
      
      hr(),
      
      # GENE PLOT
      plotOutput('plot_cpg', height = 'auto'),
      
      hr(),
      
      # INDIVIDUAL CPG PLOT
      plotOutput('boxplot_individual')
      
      #boxplot that shows tissue and trimester
      #selectInput('num','Number of CpGs to plot', choices = c(25, 50, 75,100)),
      #plotOutput('boxplot_total',  height = '650px'),
    )
  ),
  
  #footer
  hr(),
  tags$h6(paste("[1] When selecting a gene, all CpGs that are mapped to any component of an associated transcript will be shown. Transcript components include 1-5Kb upstream of the TSS, the promoter (< 1Kb upstream of the TSS), 5’UTR, exons, introns, and 3’UTRs. Everything else are intergenic regions")),
  
  tags$h6(paste("[2] All DNA methylation data here is 850k array data. See primary publication for details on processing.")),
  
  tags$h6(paste("[3] Annotations are in hg19 coordinates.")),
  
  tags$h6(paste("[4] Differential methylation is defined as a bonferroni-adjusted p-value < 0.01, and an absolute difference in mean DNAm > 25%.")),
  
  tags$h6(paste("[5] PMD regions coordinates are taken from D. I. Schroeder et al. 2013."))
  
)

server <- function(input, output, session) {
  #searching list
  updateSelectizeInput(session = session, inputId = 'name', 
                       choices = cpg$value, server = TRUE)
  
  updateSelectizeInput(session = session, inputId = 'gene', 
                       choices = gene$gene, selected = '',server = TRUE)
  
  #produce output when click submit button
  event_submit <- eventReactive(input$submit, {input$name})
  event_submit2 <- eventReactive(input$submit, {input$gene})
  empty_df <- anno[0,] %>% rename_columns() %>% select(CpG:`Relation to Transcript Features`)
  
  #cpg site information datatable
  output$info <- DT::renderDataTable({ 
    if (input$submit == FALSE)
      return(DT:: datatable(empty_df, 
                            options = list(paging = FALSE, 
                                           searching = FALSE)))
    DT::datatable(cpg_info_total(cpg_name = event_submit(), 
                                 gene_name = event_submit2())[,input$checkbox, 
                                                              drop = FALSE], 
                  options = list(pageLength = 5),
                  selection = 'single',
                  rownames = FALSE) %>% 
      formatStyle(1, cursor = 'pointer', color = 'blue')
  }
  )
  
  # find the number of transcripts for selected gene
  n_trans <- function(gene_name = NULL){
    anno_gene <- anno %>%
      filter(grepl(paste0('\\<', gene_name, '\\>'), genes_symbol),
             cpg %in% rownames(betas)) %>%
      arrange(desc(start))
    n_trans <- anno_gene$genes_tx_id %>% 
      str_split(', ') %>% 
      unlist() %>% 
      unique %>% 
      length()
    n_trans
  }
  #plot height changes depends on the number of transcripts
  plotHeight <- reactive({
    if (n_trans(event_submit2()) > 20)
    {
      return(1000)
    }
    else
    {
      return(900)
    }
  })
  #gene cpg plots
  output$plot_cpg <- renderPlot(
    {
      stat_plot(cpg_name = event_submit(), gene_name = event_submit2()) 
      
    },
    #set plot height depends on the number of transcipts
    height = plotHeight
    
  )
  
  #boxplot
  output$boxplot_total <- renderPlot(
    { 
      boxplot_total(cpg_name = event_submit(), 
                    gene_name = event_submit2(), 
                    nrow = (as.numeric(input$num) %/% 7))
    }
  )
  
  #related cpg site
  b <- reactive(fn_cpg(event_submit2()))
  #searching list (cpg sites related to selected gene)
  observe(
    {
      updateSelectizeInput(session = session, 
                           inputId = 'cpg', 
                           choices = b(), 
                           server = TRUE, 
                           selected = '')
      
    })
  
  ##retrive cpg name when click the datatable
  event_submit3 <- eventReactive(input$submit2, {input$cpg})
  
  #produce boxplots of selected cpgs
  observeEvent(input$info_cell_clicked, 
               {
                 c_info = input$info_cell_clicked
                 
                 output$boxplot_individual <- renderPlot(
                   {
                     if (!is.null(c_info$value))
                     {
                       boxplot_individual(c_info$value)
                     }
                     else if(!is.null(event_submit3()))
                     {
                       boxplot_individual(event_submit3())
                     }
                   }
                 )
               })
  
  observeEvent(input$submit2, 
               {
                 output$info <- 
                   DT::renderDataTable(
                     DT::datatable(
                       cpg_info_total(
                         cpg_name = event_submit(), 
                         gene_name = event_submit2())[,input$checkbox, drop = FALSE], 
                       options = list(pageLength = 5),
                       selection = 'single',
                       rownames = FALSE) %>% 
                       formatStyle(1, cursor = 'pointer', color = 'blue'))
               })
  
  #reset button
  observeEvent(input$reset, {reset('name')})
  observeEvent(input$reset, {reset('gene')})
  observeEvent(input$reset2,{reset('cpg')})
}
shinyApp(ui, server)

#' 
