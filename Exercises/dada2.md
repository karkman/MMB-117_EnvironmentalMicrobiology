MMB-117
================

# DADA2 pipeline

For this pipeline you need to open a RStudio session and allocate 8h,
50Gbb and 4 cores.

We will use the DADA2 to process the 16S rRNA data.  
The pipeline consists of the following steps:  
1. Quality filtering and trimming  
2. Learning the error rates  
3. Denoising the data  
4. Merging the forward and reverse reads  
5. Removing chimeras  
6. Assigning taxonomy

But before we wtart we need some setting up.

We need some R packages to run the analyses.  
And we need to set the working directory to our course folder. So change
the path to your own working directory below. And then run the code
block below.

``` r
workdir <- "PATH_TO_WORKDIR"
setwd(workdir)

.libPaths(c("/projappl/project_2013123/project_rpackages_r421", .libPaths()))
libpath <- .libPaths()[1]

library(dada2)
library(phyloseq)
library(microViz)
library(tidyverse)
```

Then we will make two lists of the read files in the folder, one for R1
and one for R2 read files. And get the samples names from the names of
the read files (from R1 read files to be precise.)

``` r
path <- paste0(workdir, "/02_TRIMMED")

fnFs <- sort(list.files(path, pattern = "_trimmed_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_trimmed_2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

To be able to decide the trimming parameters, we’ll plot the quality
profiles from first six read files. Separately for R1 and R2 reads.

``` r
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])
```

Before stating to trim the reads, we define the output files for the
trimmed reads.

``` r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Then we can quality trim the reads, but first we need to decide some the
trimming parameters based on the quality plots. Before running the
trimming command, set the trimming length (`truncLen`) and the maximum
expected errors (`maxEE`) per read.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
    truncLen = c(XXX, YYY),
    maxN = 0, maxEE = c(X, X), truncQ = 2, rm.phix = TRUE,
    compress = TRUE, multithread = TRUE
)
out
```

Next we will learn the error rates from the data and denoise the data.  
We’ll use a bit less data for this step to speed up the process (default
is 1e8).

``` r
errF <- learnErrors(filtFs, multithread = TRUE, nbases = 1e7)
errR <- learnErrors(filtRs, multithread = TRUE, nbases = 1e7)
```

And plot the error profiles.

``` r
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
```

To speed up, we’ll first dereplicate the data.

``` r
derepF <- derepFastq(filtFs, verbose = TRUE)
derepR <- derepFastq(filtRs, verbose = TRUE)
```

Then we can denoise the data. This will take some time.

``` r
dadaFs <- dada(derepF, err = errF, multithread = TRUE)
dadaRs <- dada(derepR, err = errR, multithread = TRUE)
```

And merge the forward and reverse reads.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
```

And convert it to a sequence table and print the dimensions and the
length distribution of the ASVs.

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
```

Now we can remove the chimeras from the data and check the proportion of
the data that is left.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)
```

Last part is to track the data through the pipeline.

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
```

All the ASVs are now made and we can save the sequence table as an RDS
file.

``` r
saveRDS(seqtab.nochim, "seqtab.rds")
```

Then we can assign taxonomy to the ASVs. And save the taxonomy table as
an RDS file. This will again take some time.

``` r
taxa <- assignTaxonomy(seqtab.nochim, "/scratch/project_2013123/DB/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread = TRUE)
saveRDS(taxa, "taxa.rds")
```

The DADA2 pipleine is done and we have the ASV table and the taxonomy
table saved as RDS files. We can move on to some preliminary data
visualization and exploration.

First we make a phyloseq object from the ASV and taxonomy tables. You
don’t to read the RDS files again, if you have not closed the R session.
We will also add the metadata to the phyloseq object and store the ASV
sequences as a separate data in the phyloseq object.  
The ASVs are named as their sequences at this stage and we will rename
them with ASV and a running number.

``` r
workdir <- "PATH_TO_WORKDIR"

ASV_table <- readRDS("seqtab.rds")
TAX_table <- readRDS("taxa.rds")

physeq <- phyloseq(otu_table(ASV_table, taxa_are_rows = FALSE), tax_table(as.matrix(TAX_table)))

dna <- Biostrings::DNAStringSet(taxa_names(physeq))
names(dna) <- taxa_names(physeq)
physeq <- merge_phyloseq(physeq, dna)

taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))

MMB117metadata <- read_delim(paste0(workdir, "doc/metadata.txt"), delim = "\t") %>% column_to_rownames("Sample")
sample_data(physeq) <- MMB117metadata %>% sample_data()

saveRDS(physeq, "physeq.rds")
```

Preliminary community data exploration using microViz package. Again,
you don’t need to read the RDS file again, if you have not closed the R
session.

``` r
physeq <- readRDS("physeq.rds")

ord_explore(physeq)
```
