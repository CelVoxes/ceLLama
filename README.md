ceLLama
================

![](ceLLama_files/cellama.png) ceLLama is a simple automation pipeline
for cell type annotations using large-language models (LLMs).

It has several advantages:

- Works locally, thus no information leak.
- Takes negative genes into account.
- Quite fast.
- Trained on almost all internet!

ceLLama can be super useful for quick & dirty cell type checks!

## How to install

``` r
devtools::install_github("eonurk/ceLLama")
```

## How to use

First, you need to download [`ollama`](https://ollama.com/).

Then you can choose the model of your choice. Currently, one of the best
open source LLM models is Llama3. You can run it on your terminal simply
using:

``` bash
ollama run llama3
```

This starts a local server on your machine, and you can see if it is
running by checking <http://localhost:11434/>. It should say “Ollama is
running”.

Then you are ready to go!

``` r
library(Seurat)
library(tidyverse)
library(httr)

pbmc.data <- Read10X("../../Downloads/filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# note that you can chain multiple commands together with %>%
pbmc <- SCTransform(pbmc, verbose = F) %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(dims = 1:10, verbose = F) %>%
    FindClusters(resolution = 0.5, verbose = F) %>% 
    RunUMAP(dims = 1:10, verbose = F)

DimPlot(pbmc, label = T, label.size = 3) + theme_linedraw() + theme(aspect.ratio = 1)
```

![](README_files/figure-gfm/pbmc2700-1.png)<!-- -->

``` r
# Find cluster markers
pbmc.markers <- FindAllMarkers(pbmc, verbose = F)

# split into a lists per cluster
pbmc.markers.list <- split(pbmc.markers, pbmc.markers$cluster)
```

``` r
# run cellama!
# set seed, make temperature 0 for reproducible results
library(ceLLama)

res <- ceLLama(pbmc.markers.list, temperature = 0, seed = 101)
```

    ## >> Response: Dendritic cells (DCs)

    ## >> Response: Neutrophils

    ## >> Response: CD8+ T cells

    ## >> Response: Plasma B cells

    ## >> Response: CD8+ T cells (Tc1)

    ## >> Response: NK cells

    ## >> Response: Macrophage/Monocyte

    ## >> Response: CD8+ T cells (cytotoxic/suppressor)

    ## >> Response: Macrophage/Monocyte

    ## >> Response: Myeloid cells (e.g., neutrophils or monocytes)

``` r
# transfer the labels
annotations <- map_chr(res, 1)

names(annotations) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, annotations)

DimPlot(pbmc, label = T, repel = T, label.size = 3) + theme_linedraw() + theme(aspect.ratio = 1)
```

![](README_files/figure-gfm/transfer%20annotations-1.png)<!-- -->

## Creating Reports

You can also create custom reports explaining why the annotations were
assigned.

``` r
# Get the reason for the annotation! (a bit slower)
res <- ceLLama(pbmc.markers.list, temperature = 0, seed = 101, get_reason = T)

# These creates 
generate_report_md(res)
create_html_report()
```

![](ceLLama_files/report_example.png)

You could check the full report [here](report.html).

#### Disclaimer

> LLMs make mistakes, please check important info.

## License

CC BY-NC 4.0

Please refer to <https://creativecommons.org/licenses/by-nc/4.0/>.
