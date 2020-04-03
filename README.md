## D-cellerate

#### Introduction

D-cellerate is a web application that provides a graphical interface for Seurat version 2, a popular single-cell RNA-seq (scRNA-seq) package for R. It provides an easy-to-use UI and the ability to export results as an HTML file.

#### Installing

D-cellerate is distributed as an R package, which can be installed by using the devtools package.

` devtools::install_github("https://github.com/BioData-PT/D-cellerate") `

This will install the D-cellerate package and package dependencies. 

Alternatively, you can obtain the D-cellerate Docker image by typing the following in your terminal:

`docker pull hmcs/d-cellerate-docker:latest`

The Docker image is useful to deploy D-cellerate quickly and without having to handle dependency issues.

### Launching

To launch the D-cellerate package, use the command `Dcellerate::launchApp()` in your R console. This will launch an instance of D-cellerate on port 4800. To access D-cellerate, type `localhost:4800` or `127.1.0.0.1:4800` in your web browser's address bar.

If you are using the Docker version, enter the command `docker run -d -p 3838:3838 hmcs/d-cellerate-docker` on your terminal. This command should run an instance of D-cellerate on localhost:3838

*Note for Docker Toolbox users*: If you are using the Docker Toolbox for Windows and Mac, D-cellerate will not run on localhost. Instead, it will run on port 3838 of the Toolbox virtual machine's IP, which can be found with the command docker-machine ip or in the first message that appears when booting Docker Toolbox. The default IP is usually 192.168.99.100, meaning the app will run on 192.168.99.100:3838

#### Features

D-cellerate implements a variety of features from the Seurat package with an easy to use GUI, as well as providing methods to filter data.

#### Uploading Data

D-cellerate accepts three types of data as an input.

 1. Tabular data, such as .csv or .tsv files
 2. Drop-seq tools outputs
 3. Cellranger HDF5 output files

#### Filtering Options

D-cellerate includes several options to filter your dataset.

![](https://i.gyazo.com/4c67917856591411420342d85ec5b3c9.png)

*Part of the filtering menu*

Filtering options are as follows:

- **Min UMI**: Removes cells where the number of unique molecular identifiers (UMI) found is lower than the specified number. 
- **Max UMI**: Removes cells where the number of UMI found is higher than the specified number.
- **Min Genes**: Removes cells where the number of genes found is lower than the specified number.
- **Max Genes**: Removes cells where the number of genes found is higher than the specified number.
- **Min Cells**: Removes genes based on number of cells where it is expressed.
- **Pattern for mitochondrial genes**: Subsets mitochondrial genes based on a regular expression  (e.g. ^mt- for mouse or ^MT- for human).
- **Sample cells**: Randomly samples a subset of cells that pass the above filters, allowing for quicker data exploration before using the full dataset.

Filtering results are displayed in adjacent plots, as well as a table that provides a summary of the filtering process.

The barcode plot shows the total UMI per barcode, where the barcodes that were not filtered out and were selected for further analysis are displayed in the green line.

![](https://i.gyazo.com/917898911ad0151da7fc169b41556804.png)

*Barcode plot*

A set of violin plots shows the total distributions of total UMI per cell (left plot), and the number of detected genes per cell (right plot).

![](https://i.gyazo.com/76b298446f493af85e58ed725587c4c0.png)

*Violin plots*

A summary allows the user to track the filtering process, and see how many barcodes and genes passed through for further analysis.

![](https://i.gyazo.com/131aaca882fe9fed5539cf3652141d8c.png)

*Summary*

#### Normalization

Normalizes gene expression across all cells through log transformation. The output is presented in a barplot.

![](https://i.gyazo.com/9ae9bc93f31860c537f514cf34cbf50b.png)

*Log-normalization menu*

![](https://i.gyazo.com/c1da8a013497906a092f439e13fac052.png)

*Histogram of mean gene expressin across all cells*


#### Variable Genes

D-cellerate allows the user to identify variable and non-variable genes, as well as identify features that are outliers on a mean variability plot.

![](https://i.gyazo.com/655a208be245a4823ccb7f1a94541f34.png)

*Variable genes settings*

![](https://i.gyazo.com/fc868c1c10f521efdd548cbc150d91be.png)

*Mean variablity plot*


#### Dimensional Reduction and Clustering

D-cellerate performs a principal components analysis (PCA) with either the variable genes or all genes in the cell. It will output a scree plot showing the proportion of variance explained by each principal component, as well as Scatterplots of the output of the principal component analysis. 

![](https://i.gyazo.com/40ca517d42491707829f576053a7e1da.png)

*PCA plot*

![](https://i.gyazo.com/82226dd8d0b308c7482531073e6cb26e.png)

*PCA options menu*

![](https://i.gyazo.com/48eb17fc88d9da33baa21a3bdaf27de5.png)

*PCA plots*

![](https://i.gyazo.com/38525e0e5cd28efcf9c69ab8cc4b9d0c.png)

*Elbow plot*

![](https://i.gyazo.com/441a46f04e4c277d24953cc70b436a9a.png)

Dimensional reduction is performed using t-SNE, and is followed by a clustering analysis.

![](https://i.gyazo.com/40d48ada220b5996d5ded76687758e4f.png)

*t-SNE options*

![](https://i.gyazo.com/6cb49dea2f7f9037fa3f869adb5a3083.png)

*Cluster plot*

#### Differential Expression and Marker Gene Identification

Marker gene discovery tests the various clusters generated in the clustering test to find markers. The available tests are the Wilcoxon rank sum, Student's t-test and Student's AUC classifier.

![](https://i.gyazo.com/eea711c6c58552e8b0952cde836fb60a.png)

*Differential expression menu*

There are three outputs for this section: A differential expression table, a gene expression heatmap and t-SNE cluster plots for marker visualization in clusters.

The differential expression table finds differentially expressed genes between clusters, and displays them in a tabular format.

![](https://i.gyazo.com/8be5b503f86ee218ea3c5be242385cf1.png)

*Differential expression table*

The heatmap shows the top markers for each cluster.

![](https://i.gyazo.com/a5f2cae3dfb8cc005236c61a4bc60e75.png)

*Top markers heatmap*

The t-SNE plot(s) show marker gene expression per cluster.

![](https://i.gyazo.com/619179899bee8256d36c00b66d26e977.png)

*Marker visualization plots*

#### Export Analysis

D-cellerate exports results as an HTML report. In order to ensure analysis reproducibility, the report includes both plots and code. Furthermore, if the user intends to perform various analyses, the user can bookmark the application state, saving the variables used for the analysis so that another analysis can be safely run.
