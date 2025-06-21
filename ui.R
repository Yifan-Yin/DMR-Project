
#UI ----
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
      radioButtons('checkbox_trimester',
                   'Include samples (trimester):',
                   setNames(c('First', 'Third'), c('First Trimester', 'Third Trimester (Term)')),
                   selected = 'Third'),
      checkboxGroupInput('checkbox_tissue',
                         'Include samples (cell type):',
                         setNames(unique(pDat$Tissue), gsub(' cs', '', unique(pDat$Tissue))),
                         selected = unique(pDat$Tissue)),
      
      
      #checkboxGroupInput('checkbox_sex',
      #                   'Include samples (sex):',
      #                   setNames(c('F', 'M'), c('Female', 'Male')),
      #                   selected = unique(pDat$Sex)),
      
      
      
      sliderInput('range_cpg', 
                  'Select number of CpGs to plot (only applies to Condensed View):',
                  min = 1, max = 500, value = 500),
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
                             'DNA methylation',
                             min = 0, max = 25, value = c(10), step = 1)
                         ),
                         column(
                           3,
                           sliderInput(
                             'condensed_plot_height_2',
                             'Difference from CpG mean',
                             min = 0, max = 25, value = c(5), step = 1)),
                         column(
                           3, 
                           sliderInput(
                             'condensed_plot_height_3',
                             'Annotations',
                             min = 0, max = 25, value = c(5), step = 1)),
                         column(
                           3, 
                           sliderInput(
                             'condensed_plot_height_4',
                             'Transcripts',
                             min = 0, max = 25, value = c(5), step = 1)
                         )
                       )
                     ))
                 ),
                 
                 br(),
                 textOutput('text'),
                 plotOutput('annoplot', height = 900)),
        
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
                             'expanded_plot_height_1', 'DNA methylation',
                             min = 1, max = 10, value = c(3), step = 1)
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
                 plotOutput('expanded', height = 900))
        
        #tabPanel("CpG Boxplots",
        #        p('Show boxplots for specific CpGs:'),
        
        #         cpgs search list related to the input
        #         selectizeInput('cpg', 'CpGs related to the input',
        #                        choices = NULL,
        #                        multiple = TRUE),
        
        #         action button to submit and reset input
        #         actionButton('submit_cpg_gene_input', 'Submit'),
        #         actionButton('reset_cpg_gene_input','Reset'),
        
        #         plotOutput('boxplot_selected_cpg'))
        
      ))),
  
  tags$h3('Annotation Table'),
  fluidRow(
    div(
      dataTableOutput('dt_anno_cpg', width = "98%"), 
      style = "font-size:85%",
      align = 'center')),
  
  #footer
  hr(),
  
  HTML(paste(
    "[1] Only associated CpGs will be displayed when selecting a gene.This CpG-gene mapping is based on the Illumina-provided EPIC array annotation.",
    "[2] When selecting a gene, all CpGs that are mapped to any component of an associated transcript will be shown. Transcript components include 1-5Kb upstream of the TSS, the promoter (< 1Kb upstream of the TSS), 5’UTR, exons, introns, and 3’UTRs. Everything else are intergenic regions", 
    "[3] All DNA methylation data here is 850k array data. See primary publication for details on processing.",
    "[4] Annotations are in hg19 coordinates.",
    "[5] Differential methylation is defined as a bonferroni-adjusted p-value < 0.01, and an absolute difference in mean DNAm > 25%.",
    "[6] PMD regions coordinates are taken from D. I. Schroeder et al. 2013.",
    "[7] Differential methylation by cell type was conducted separated in first and third trimester samples for autosomal CpGs. Differentially methylated cytosines (DMCs) are defined as a bonferroni adjusted p-value < 0.01 and a mean difference in methylation > 10%. Differential methylation was conducted with a one-versus-all design.",
    sep = " <br> "))
)
