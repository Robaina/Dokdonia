# :bulb: About

This repo contains the source code employed in the study "Gene expression response to temperature acclimation in a model photoheterotrophic marine flavobacterium". Source code is organized into a Python package, and dependencies are managed with Conda.

## :wrench: Installation

To install the package, clone the repo and run the setup script:

```bash
git clone https://github.com/Robaina/Dokdonia.git
cd Dokdonia
bash setup.sh
```

## :rocket: Usage

To use the package, activate the Conda environment:

```bash
conda activate dokdonia
```

This environment includes a Jupyter Notebook ipykernel, so the package can be used in a notebook as well.

## :notebook_with_decorative_cover: Notebooks

The steps followed to generate the results and figures in the paper are detailed in the following notebooks:

- [Differential expression](notebooks/1_differential_expression.ipynb)
- [Clustering: DeSeq2-normalized counts](notebooks/2_clustering_deseq2.ipynb)
- [Clustering: Transcript abundances](notebooks/3_clustering_tc.ipynb)
- [Figures](notebooks/4_paper_figures.ipynb)