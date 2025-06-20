DMR-Project
================

A *brief* introduction of what is in this repo…

| File                  | Description                              |
| :-------------------- | ---------------------------------------- |
| shiny\_DMR.Rmd        | code to run app                                         |
| updated\_plots.Rmd    | Code for updated plots. October 11, 2019 |


The app depends on these inputs:

# sample data
pDat <- readRDS(here('data-production', 'pDat.rds'))

# methylation data
betas <- readRDS(here('data-production', 'betas_summarized.rds'))

# annotation
anno <- readRDS(here('data-production', 'annotation.rds'))

#color code
color_code <- readRDS(here::here('data-production', '2_3_color_code.rds'))
color_code_tissue <- setNames(color_code$Colors_Tissue, color_code$label)
colors <- color_code_tissue[unique(pDat$Tissue)]

#  transcript models for exapnded view
anno_annotatr <- readRDS(here('data-production', 'annotation_annotatr.rds'))