# Guided selection of single cell clustering parameters through sub-sampling cluster robustness metrics

This repository contains an example implementation in `R` using `Seurat` of the framework outlined in: 

> Patterson-Cross, R.B., Levin, A.J. & Menon, V. Guided selection of single cell clustering parameters through sub-sampling cluster robustness metrics. (2020).

## Installation

To use the code and examples in the repository, first clone the repository to your computer:

```bash
git clone https://github.com/rbpatt2019/cluster.stability.git
```

The code is implemeted in R, and the dependencies are pinned in the `renv.lock` file. To install dependencies, open an R terminal, then proceed as follows:

```R
install.packages("renv") # Not necessary if already installed
renv::init()
```

This will create a local project and install the dependencies there, rather than into your global `R` installation.

## Usage

Two example scripts are included in this repo. The [first](./examples/1_seurat_pipeline.R) runs through the main analysis framework and covers the key steps, including iterative, sub-sampled clustering, calculating the co-clustering frequency matrix, determining the silhouette scores, and creating the silhouette distribution plots.
The [second](./examples/2_seurat_further_visuaisations.R) covers several additionaly visualisations that we used in the paper and find useful for understanding the patterns and clusters within your data. 
Both these files can be run as stand-alone scripts, as so:

```R
Rscript examples/1_seurat_pipeline.R
Rscript examples/2_seurat_further_visuaisations.R
```

If used this way, they are meant to be run sequentially. Alternatively, they can be opened with any modern IDE or editor for interactive execution.

Additionally, the interested user may directly source the scripts in the `R` directory to create analyses that suit individual needs.

## Data

An example data set is included with this repo. It contains ~500 human PBMCs sequenced with Smart-Seq from [Ding, et al. _Nature Biotechnology_, 2020](https://www.nature.com/articles/s41587-020-0465-8) and is one the datasets detailed in our paper. The RDS included in the `data` directory has already been normalised and had PCAs calculated using Seurat's `SCTransform` and `RunPCA` functions. If you are using your own data, be sure to normalise and calculate PCAs first. 
