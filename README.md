Placental Cell Methylome Browser
================

Current URL (as of 6/23/2025): 

- https://wvictor.shinyapps.io/dmr-project/ 

Prior URLs:

- (pre 2025): https://robinsonlab.shinyapps.io/Placental_Methylome_Browser/



The app depends on these inputs:

#### sample data
`pDat <- readRDS(here('data-production', 'pDat.rds'))`

`pDat` contains the sample metadata

#### methylation data

`betas <- readRDS(here('data-production', 'betas_summarized.rds'))`

#### annotation

Victor's custom annotation for CpGs

`anno <- readRDS(here('data-production', 'annotation.rds'))`

#### color code

For visualizations

```
color_code <- readRDS(here::here('data-production', '2_3_color_code.rds'))
color_code_tissue <- setNames(color_code$Colors_Tissue, color_code$label)
colors <- color_code_tissue[unique(pDat$Tissue)]
```

####  transcript models for exapnded view

For expanded view need the overlapping transcript models 

`anno_annotatr <- readRDS(here('data-production', 'annotation_annotatr.rds'))`
