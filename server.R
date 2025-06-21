server <- function(input, output, session) {
  ############  searching lists ############
  #searching list of cpgs
  updateSelectizeInput(session = session, inputId = 'cpg_name', 
                       choices = cpg$value, server = TRUE)
  #searching list of genes
  updateSelectizeInput(session = session, inputId = 'gene_name', 
                       choices = gene$gene, selected = '',server = TRUE)
  #searching list of chromosomes
  chr <- unique(anno$chr) %>% str_sort(numeric = TRUE) %>% as_tibble()
  updateSelectizeInput(session = session, inputId = 'chr',
                       choices = chr$value, selected = '', server = TRUE)
  
  #get the inputs and store as reactive values 
  get_gene_input <- eventReactive(input$submit_gene_input, {input$gene_name})
  get_cpg_input <- eventReactive(input$submit_cpg_input, {input$cpg_name})
  get_region_chr<- eventReactive(input$submit_region_input, {input$chr})
  get_region_start<- eventReactive(input$submit_region_input, {input$start})
  get_region_end<- eventReactive(input$submit_region_input, {input$end})
  
  #get cpg sites related to the selected gene and store as reactive values
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
        make_datatable(gene_name = get_gene_input(), 
                       cpg_name = NULL,
                       Chr = NULL, Start = NULL, End = NULL)[,input$checkbox,drop = FALSE],
        options = list(pageLength = 5),
        selection = 'single',
        rownames = FALSE) %>%
        formatStyle(1,color = 'blue')}
      
      #output a datatable when input cpgs
      else if(input$tabs == 'CpG Input') {
        
        DT::datatable(
          make_datatable(gene_name = NULL, cpg_name = input$cpg_name, 
                         Chr = NULL, Start = NULL, End = NULL)[,input$checkbox, drop = FALSE], 
          options = list(pageLength = 5),
          selection = 'single',
          rownames = FALSE) %>%
          formatStyle(1, color = 'blue')
      } else if(input$tabs == 'Region Input') {#output a datatable when input a region
        
        DT::datatable(
          make_datatable(
            gene_name = NULL, cpg_name = NULL,
            Chr = get_region_chr(), Start = get_region_start(), End = get_region_end()
          )[,input$checkbox, drop = FALSE],
          options = list(pageLength = 5),
          selection = 'single',
          rownames = FALSE
        ) %>%
          formatStyle(1, color = 'blue')
        
      })
  }
  )
  
  ############  annoplot output: default to show first 250 cpgs############
  #annotation_track plots
  observeEvent(
    c(input$submit_region_input,
      input$submit_gene_input),
    {
      
      output$annoplot <-renderPlot(
        {
          #return empty when no inputs
          if (input_value$input_value  == FALSE) return()
          #output a plot when input a gene
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
                                 Start = get_region_start(), End = get_region_end(), 
                                 Chr = get_region_chr(), 
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
    }
  )
  ############  annoplot output: choose cpg number############
  #store cpg number into reactive value
  cpg_number <- reactiveValues(cpg_number = 250)
  #get the cpg number input
  observeEvent(input$submit_cpg_set, {cpg_number$cpg_number <- input$range_cpg})
  #reset the number of cpgs to plot into default(250) when update inputs
  observeEvent(input$submit_gene_input, {cpg_number$cpg_number = 250})
  observeEvent(input$submit_region_input, {cpg_number$cpg_number = 250})
  #update the plot when change the cpg number to show
  observeEvent(input$submit_cpg_set, {output$annoplot <- renderPlot(
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
                           cpg_number = cpg_number$cpg_number, 
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
        expanded_view(gene_name = get_gene_input(), 
                      
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
  #get_cpg_from_dt <- eventReactive(input$submit_cpg_gene_input, {input$cpg})
  
  #produce a boxplot based on the selected cpg from datatable
  #observeEvent(input$dt_anno_cpg_cell_clicked, 
  #             {
  #               c_info = input$dt_anno_cpg_cell_clicked
  
  #               output$boxplot_selected_cpg <- renderPlot(
  #  {
  #    if (!is.null(c_info$value))
  #    {
  #      boxplot_selected(c_info$value, sex_type = input$checkbox_sex)
  #    }
  #    else if(!is.null(get_cpg_from_dt()))
  #    {
  #      boxplot_selected(get_cpg_from_dt(), sex_type = input$checkbox_sex)
  #    }
  #  }
  #)
  #})
  #produce a boxplot based on the selected cpg using search box
  #observeEvent(input$submit_cpg_gene_input, {output$boxplot_selected_cpg <- renderPlot(
  #  {
  #    boxplot_selected(get_cpg_from_dt(), sex_type = input$checkbox_sex)
  #  }
  #)})
  
  ############  reset boxplot output ############
  #reset the associated cpg inputs in search box
  #observeEvent(input$reset_cpg_gene_input, {reset('cpg')})
  #reset the boxplot when click the reset button
  #observeEvent(input$reset_cpg_gene_input, {output$boxplot_selected_cpg <- renderPlot(
  #  {
  #    return()
  #  }
  #)})
  #reset the boxplot when switch between tabs
  #observeEvent(input$tabs, {output$boxplot_selected_cpg <- renderPlot(
  #  {
  #    return()
  #  }
  #)})
  #reset the boxplot when change into another gene
  #observeEvent(input$submit_gene_input, {output$boxplot_selected_cpg <- renderPlot(
  #  {
  #    return()
  #  }
  #)})
  #reset the boxplot when reset gene input
  #observeEvent(input$reset_gene_input, {output$boxplot_selected_cpg <- renderPlot(
  #  {
  #    return()
  #  }
  #)})
  
  
}