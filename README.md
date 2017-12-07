# GENIE3

## Description
GENIE3 (GEne Network Inference with Ensemble of trees): Inference of gene regulatory networks from gene expression data

GENIE3 is now available from [Bioconductor](https://bioconductor.org/packages/devel/bioc/html/GENIE3.html). A **tutorial** [(vignette)](https://bioconductor.org/packages/devel/bioc/vignettes/GENIE3/inst/doc/GENIE3.html) is included in the package.

## Official

### Installation
To install GENIE3 (from Bioconductor development branch), you can run the following commands from R:
```
install.packages("https://bioconductor.org/packages/devel/bioc/src/contrib/GENIE3_0.99.7.tar.gz", repos=NULL)

# You might need to install these packages first:
install.packages(c("stats", "reshape2"))
```

## This branch

Features added:

* multi-node parallelism
* progress bar

### Installation
To install GENIE3, you can run the following commands from R:
```
install_git("git://github.com/mase5/GENIE3.git", branch = "mnp_pb_features")

# You might need to install these packages first:
install.packages(c("stats", "reshape2", "doSNOW", "doRNG", "foreach", "doParallel"))
```

### How to run

Here is a snippet:
```
# Define the cluster configuration
psc<-list(config=list(type="raw", "def"=list(user = "userid", nodes = c("[machine-address-1]","[machine-address-2]"), n.cores = c([n-cores-1],[n-cores-2]), verbose = T)))
weight.matrix<-GENIE3(data.matrix, regulators=inputTFs, cluster = psc, progress = T, verbose=TRUE)
```
NOTE: [machine-address-1] should be the address where the R script is started

![alt text](https://raw.githubusercontent.com/mase5/GENIE3/mnp_pb_features/images/progress.gif)


