# go-cluster
Clustering PPINs by GO similarity scores.

Uses [**STRINGdb**](https://string-db.org/) and [**INet**](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190029) networks, [**GSE23561**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23561) microarray gene expressions data, [**GOSemSim**](https://bioconductor.org/packages/release/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html) and [**org.Hs.ed.db**](https://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf) packages.

## Running instructions:
You first need download and export link files into ``LINKS/`` directory using the download link specified in ``LINKS/DOWNLOAD.txt``.

### Dependencies:
To overcome possible dependency problems, run scripts in this order:

``degs.R`` :arrow_right: ``links.R`` :arrow_right: ``msLinks.R`` :arrow_right: ``go.R``

> **Note:** All scripts are dependent to ``vars.R``
