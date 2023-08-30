# scdney - Single cell data integrative analysis

<!--
<br />
-->


<!--
[![](https://img.shields.io/github/last-commit/SydneyBioX/scdney.svg)](https://github.com/SydneyBioX/scdney/commits/master) [![Travis build status](https://travis-ci.org/SydneyBioX/scdney.svg?branch=master)](https://travis-ci.org/SydneyBioX/scdney)


<br />

-->




`scdney` is a suite of packages single cell analysis R packages developed by members of the Bioinformatics Cluster within the [Sydney Precision Data Science Center](https://www.sydney.edu.au/science/data-science) at The University of Sydney.



<link rel="stylesheet" href="css/hexagons.css">

<div id="importedContent"></div>
    <script>
        fetch('https://raw.githubusercontent.com/SydneyBioX/packageHeader/main/spatialHeader.html')
            .then(response => response.text())
            .then(htmlContent => {
                const importedContentDiv = document.getElementById('importedContent');
                importedContentDiv.innerHTML = htmlContent;
            })
            .catch(error => {
                console.error('Error fetching content:', error);
            });
    </script>



## Installation

Before running the installation command below, some `R` dependencies must be installed. This can be done as follows:

```r
install.packages(c("BiocManager","remotes"))
```

To install all packages in available `scdney` and attached the installed `scdney` packages in the current R session.

```r
# install scdney packages
BiocManager::install("SydneyBioX/scdney")
# load  scdney packages
library(scdney)
```

### Installation Notes
#### Installation Time

Installing the whole suite of packages can take a long time on UNIX-based systems:

+ Around 6 mins when using prebuilt binaries from [r2u](https://eddelbuettel.github.io/r2u/) (discussed below). Unfortunately, this only works on Ubuntu LTS systems.
+ Around 20 mins on a basic r installation including some basic Bioconductor packages (e.g., [bioconductor_docker](https://github.com/bioconductor/bioconductor_docker) image).
+ More than 1 hour on system with a basic r installation (e.g., [r-base](https://hub.docker.com/_/r-base/) image)

If compiling the binaries locally is a priority, this issue can be partially addressed by adding `MAKE='make -j NCORES'` to one's `Renviron` file.

+ `NCORES` should be replaced be the number of jobs to run in parallel.

Alternatively, this process can be made significantly faster by using [`r2u`](https://eddelbuettel.github.io/r2u/). This project builds `r`
 packages for Ubuntu LTS repository system and integrates them with the system package manager. This will lead to significant increases in performance. An example docker build image for this method is included in the `/build` directory.
 
#### Non-R Dependencies

There are several non-R dependencies that must be installed in order for `scdney` to install correctly. Most of these are installed by default, but if there is a problem in the installation, one of this system requirements is likely to be the culprit.

 + `make`, `libpng-dev`, `pandoc`, `libjpeg-dev`, `zlib1g-dev`, `libfreetype6-dev`, `libfribidi-dev`, `libharfbuzz-dev`, `libxml2-dev`, `libfontconfig1-dev`, `pandoc-citeproc`, `cmake`, `libicu-dev`, `libssl-dev`, `libglpk-dev`, `libgmp3-dev`, `libtiff-dev`, `libcurl4-openssl-dev`



## `scdney` Workshops Series

+ [**Single-cell analysis fundation workshop**](https://sydneybiox.github.io/BIS2019_SC/), where we use two mouse liver datasets to illustrate three critical topics in scRNA-seq analysis: Quality control of scRNA-seq data, data integration of multiple scRNA-seq data,
and cell type marker identification.
+ [**Advanced Phenotyping using scdney**](https://sydneybiox.github.io/scdneyAdvancedPhenotyping/articles/advanced_phenotyping.html), where we describe various strategies for cellular phenotyping, construct and characterise trajectory, and identify cell-cell interaction.
+ [**CITE-seq data analysis using CiteFuse**](https://sydneybiox.github.io/BiocAsia2020CiteFuse/articles/CiteFuse_BioCAsia_workshop.html), a hands-on experience to the CiteFuse package workshop.
+ [**Disease Outcome Classification in Single Cell Patient Data Analysis**](https://sydneybiox.github.io/scdneyDiseasePrediction/articles/disease_outcome_classification_schulte.html), where we explore various strategies for disease outcome prediction using single cell data.

## Citation

Cao Y, Tran A, Kim H *et al.* Thinking process templates for constructing data stories with SCDNEY [version 1; peer review: 1 approved]. *F1000Research* 2023, **12**:261 (<https://doi.org/10.12688/f1000research.130623.1>)


