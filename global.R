#set up
library(rsconnect)
library(shiny)
# rsconnect::setAccountInfo(name='yifan7', token='52E1037550A7AA68252579A98EE8AA38', secret='+NlrJysSVZHa5LQQb9slmYRDXExPfjF0QbxMjUH1')
library(DT)
library(shinyjs)
library(shinythemes)
library(here)

library(tidyverse)
library(egg)
library(here)
library(patchwork)

source(here::here('R', 'visualization.R'))
source(here::here('R', 'server.R'))

# sample data
pDat <- readRDS(here('data', 'pDat.rds'))
# methylation data
betas <- readRDS(here('data', 'betas_summarized.rds'))
# annotation
anno <- readRDS(here('data', 'annotation.rds'))
#color code
color_code <- readRDS(here::here('data', '2_3_color_code.rds'))
color_code_tissue <- setNames(color_code$Colors_Tissue, color_code$label)
colors <- color_code_tissue[unique(pDat$Tissue)]

#  transcript models for exapnded view
anno_annotatr <- readRDS(here('data', 'annotation_annotatr.rds'))

#list of genes for search list
gene <- tibble(gene = str_split(anno$genes_symbol, ', ') %>%
                 unlist() %>%
                 unique() %>%
                 na.omit() %>%
                 sort())

#pull out cpg name for search list use
cpg <- unique(anno$cpg) 
cpg <- tibble(value = cpg[cpg %in% rownames(betas)])
