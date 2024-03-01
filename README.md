# go-cluster

Clustering PPINs by GO similarity scores.

Uses [**STRINGdb**](https://string-db.org/) and [**INet**](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190029) networks, [**GSE23561**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23561) microarray gene expressions data, [**GOSemSim**](https://bioconductor.org/packages/release/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html) and [**org.Hs.ed.db**](https://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf) packages.

## Running instructions
You first need download and export the PPI data (i.e. link files) into ``LINKS/`` directory using the download instructions given in ``LINKS/DOWNLOAD.txt``.

**To be able run the pipeline with different Gene Expression and PPI data, you must apply the following changes:**

1) Replace the ``RAW.csv`` file with your own gene expression data. This is a **comma-separated matrix**, where each column represents genes and each column represents the subjects (i.e., persons or test-conditions). Please, see the ``RAW.csv`` as an example.

2) Put the PPI file to the ``LINKS/`` directory. This file must consist of a **comma-separated edgelist** representing the PPI score (i.e., weight) for the gene pairs. Please, see the ``LINKS/elinks_inet.csv`` as an example.

3) Update the variables in the ``vars.R`` file to match your Gene Expression and PPI data. It is necessary to reassign the ranges (i.e., ``intervals``) in order to correctly separate the control group and each disease group in Gene Expression matrix. It is also necessary to reset the ``cutoff`` value for PPI scores used in PPIN reduction. In ``vars.R``, you can also change the t-test significance (``P_VAL``) and the fold-change threshold (``FC``) that are used in the DEG analysis.


## Pipeline

<a target="_blank" rel="noopener noreferrer" href="https://raw.githubusercontent.com/smtnkc/go-cluster/master/pipeline.png"><img src="https://raw.githubusercontent.com/smtnkc/go-cluster/master/pipeline.png" alt="Pipeline" width="500px"></a>

## Scripts

* ``vars.R`` manages the packages, global variables, paths, and I/O files.
* ``degs.R`` includes the necessary functions to identify differentially expressed genes.
* ``links.R`` handles reading, preprocessing, and mapping of PPI data.
* ``msLinks.R`` responsible from filtering out unmapped and insignificant PPIs.
* ``go.R`` fetches GO information and calculates the GO similarity scores.
* ``spici.R``, ``mcl.R``, ``linkcomm.R`` performs the clustering operations.
* ``bhi.R`` calculates Biological Homogeneity Index for the clusters.
* ``stability.R`` calculates stability of disease modules.
* ``validation.R`` searches for the identified genes in DEGS of the validation sets.
* ``plot.R`` responsible from generating the plots for the obtained results.
* ``stats.R``, ``bhi_stats.R`` creates tables including statistics about the clustering or validation.


## Dependencies:
It is highly recommended to install all packages required in the ``vars.R`` file.

To overcome possible dependency problems, please run the scripts in the following order, and note that all scripts depend to ``vars.R`` which manages the packages as well as the global paths, files, and variables:

``vars.R`` :arrow_right: ``degs.R`` :arrow_right: ``links.R`` :arrow_right: ``msLinks.R`` :arrow_right: ``go.R`` :arrow_right: ``spici.R || mcl.R || linkcomm.R`` :arrow_right: ``bhi.R || stability.R`` :arrow_right: ``validation.R`` :arrow_right: ``plot.R || stats.R || bhi_stats.R``

## Citation

S. Tenekeci, S. Tekir, Identifying promoter and enhancer sequences by graph convolutional networks, Computational Biology and Chemistry (2024) https://doi.org/10.1016/j.compbiolchem.2024.108040
