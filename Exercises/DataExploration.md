MMB-117
================

# Data Exploration

## Libraries

``` r
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggrepel)
```

## Setup

``` r
setwd("PATH/TO/COURSE/FOLDER")
```

Read in the phyloseq object and extract the metadata from the object.

``` r
physeq <- readRDS("physeq.rds")

MMB117metadata <- sample_data(physeq) %>%
    as_tibble() %>%
    mutate(Sample = rownames(sample_data(physeq))) %>%
    as.data.frame()

View(MMB117metadata)
```

Extract the ASV table and taxonomic table from the phyloseq object.

``` r
ASV_table <- physeq %>%
    otu_table() %>%
    as.data.frame()

tax_table <- physeq %>%
    tax_table() %>%
    as.data.frame()
```

Have a look at the tables.

``` r
str(ASV_table)
View(ASV_table)

str(tax_table)
View(tax_table)
```

Check the ASV names in both tables and see if they match.

``` r
tax_table %>%
    rownames() %>%
    head()

ASV_table %>%
    colnames() %>%
    head()

identical(rownames(tax_table), colnames(ASV_table))
```

Store the row sums (total ASV count) for later use and check the number
of ASVs that are detected only less than 3 times. We can also make a
barplot of the ASV counts per sample.

``` r
ASV_counts_per_sampleDF <- data.frame(LibSize=rowSums(ASV_table))
barplot(ASV_counts_per_sampleDF$LibSize, names=rownames(ASV_counts_per_sampleDF), main="ASV counts per sample", ylab="ASV count", las=2)

dim(ASV_table)

dim(ASV_table[, which(colSums(ASV_table) < 3)])
```

One could remove these
(`ASV_table_over3<-ASV_table[, which(colSums(ASV_table)>3)]`).  
Or singletons, but we will keep them for now.

# Normalization

We will normalize the ASV table by total sum scaling (TSS) also known as
relative abundances. And then we transpose the table so that the ASVs
are in rows and samples in columns.

``` r
ASV_tableRA.mat <- apply(ASV_table, 1, function(i) i / sum(i))
ASV_tableRA.mat <- t(ASV_tableRA.mat)

View(ASV_tableRA.mat)
```

Relative abundances should add up to 1. We can check this.  
Two different ways to do this that should give the same result.

``` r
apply(ASV_tableRA.mat, 1, sum)

rowSums(ASV_tableRA.mat)
```

We can remove the BI negative control as there is no data. And we should
remove it from the metadata as well.

``` r
ASV_tableRA.mat <- ASV_tableRA.mat[2:nrow(ASV_tableRA.mat), ]
MMB117metadata <- MMB117metadata[-1, ]
```

# Ordination

Ordination with all samples and without our own negative control.

``` r
plot(metaMDS(ASV_tableRA.mat), type = "text", display = "sites")
plot(metaMDS(ASV_tableRA.mat[2:25, ]), type = "text", display = "sites")
```

Looks pretty good!

## Diversity

We can calculate the Shannon diversity index and put that into metadata

``` r
MMB117metadata$ASV_divSha <- diversity(ASV_tableRA.mat, index = "shannon")
```

## Data exploration

According to: Zuur, A.F., Ieno, E.N. and Elphick, C.S. (2010), A
protocol for data exploration to avoid common statistical problems.
[Methods in Ecology and Evolution, 1:
3-14](https://doi.org/10.1111/j.2041-210X.2009.00001.x)

I want to start from ***INDEPENDENCE***, since it is the most important
assumption:

We had only one sampling point. Thus we don’t really have “traditional”
independence problem which would arise from repeated measurements.  
If we would have taken samples before and after some treatment from the
same sites, we would have to take this into account in our statistical
analyses (use mixed models, etc.)

### Outliers in X & Y

X = metadata and things that are usually on x-axis

Y = the numerical data we use, in this case we have microbiome data that
are either counts (non-normalized),  
relative abundance (TSS-normalized), or clr-transformed normally
distributed numbers that can be used in statistical analyses

Metadata first.

``` r
str(MMB117metadata)
```

So lets plot all variables that make sense to be plotted.

Snow depth.

``` r
ggplot(MMB117metadata, aes(Sample, SnowDepth, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("Snow depth")
```

Not recorded for all, maybe can be left out from the analysis

Wet weight.

``` r
ggplot(MMB117metadata, aes(Sample, WetWeight, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("Wet weight (g)")
```

Dry weight.

``` r
ggplot(MMB117metadata, aes(Sample, DryWeight, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("Dry weight (air dried) (g)")
```

80C dry weight.

``` r
ggplot(MMB117metadata, aes(Sample, X80C_DryWeight, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("Owen dry weight, pore water dried out (g)")
```

DNA weight.

``` r
ggplot(MMB117metadata, aes(Sample, DNAWeight, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("What was weighed for DNA extraction (mg)")
```

Moisture.

``` r
ggplot(MMB117metadata, aes(Sample, Moisture, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("Soil water % (?)")
```

Gravimetric water content would be better than “moisture”, you could
change the column. There are some that are very wet. Interesting to see
if their microbiome is completely different because of the snow for
instance

pH’s.

``` r
ggplot(MMB117metadata, aes(Sample, pH_W, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("Water-extraxted pH")

ggplot(MMB117metadata, aes(Sample, pH_Ca, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("CaCl2-extraxted pH")
```

Nice grouping!!  
With both pH’s actually

SOM.

``` r
ggplot(MMB117metadata, aes(Sample, SOM, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("Soil organic matter (LOI)")
```

Few samples that are under or above but not too extreme I would say Some
missing data as well

Community diversity.

``` r
ggplot(MMB117metadata, aes(Sample, ASV_divSha, label = Site)) +
    geom_point(color = "black", size = 3) +
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 270, hjust = 0, vjust = 0.5)) +
    geom_text_repel(color = "black", size = 4, max.overlaps = Inf) +
    ggtitle("Shannon's diversity index, Y data")
```

Negative control is an outlier, as it should

The NMDS plot also shows if there are otliers in the samples (in terms
of Y data)

``` r
plot(metaMDS(ASV_tableRA.mat[2:25, ]), type = "text", display = "sites")
```

And not really

``` r
plot(metaMDS(ASV_tableRA.mat[2:25, ]))
```

This shows both the species and sites and from this we can see that
different species associate to different sites

**Good! looks like there are no weird extreme cases in the data.**

### Homogenity of variances

With statistical analyses we usually compare means, so if the black
lines are in different levels, that is fine.  
But the variability of the observations should be similar. i.e. the size
of the boxes should be about the same

So lets check this with metadata and Shannons diversity index

``` r
ggplot(subset(MMB117metadata, Sample != "Neg. control"), aes(Site, Moisture)) +
    geom_boxplot() +
    theme_bw(base_size = 14) +
    theme(panel.border = element_rect(colour = "black")) +
    theme(legend.position = "none")
```

whoops! Good thing we are not studying change of soil water % in
response to land use! But I wouldn’t build a linear model using soil
water % as an explaining variable. Non-parametric statistical test would
be probably ok.

``` r
ggplot(subset(MMB117metadata, Sample != "Neg. control"), aes(Site, pH_W)) +
    geom_boxplot() +
    theme_bw(base_size = 14) +
    theme(panel.border = element_rect(colour = "black")) +
    theme(legend.position = "none")
```

Same here but less extreme.

``` r
ggplot(subset(MMB117metadata, Sample != "Neg. control"), aes(Site, pH_Ca)) +
    geom_boxplot() +
    theme_bw(base_size = 14) +
    theme(panel.border = element_rect(colour = "black")) +
    theme(legend.position = "none")
```

This could work with liner model as an X-variable

``` r
ggplot(subset(MMB117metadata, Sample != "Neg. control"), aes(Site, SOM)) +
    geom_boxplot() +
    theme_bw(base_size = 14) +
    theme(panel.border = element_rect(colour = "black")) +
    theme(legend.position = "none")
```

Some samples don’t have SOM -values.  
Gas station is “in another planet”, so SOM as an X-variable of a linear
model would probably be a bad idea.

``` r
ggplot(subset(MMB117metadata, Sample != "Neg. control"), aes(Site, ASV_divSha)) +
    geom_boxplot() +
    theme_bw(base_size = 14) +
    theme(panel.border = element_rect(colour = "black")) +
    theme(legend.position = "none")
```

Also diversity shows heterogeneity in the variances, but this probably
would work in linear models as an Y-variable.  
We can take these into account by selecting statistical models that that
do not require homogeneity.  
If we would really need to use linear models, we could also transform
these variables or delete observations that are outliers and go with for
instance 3 observations per site.  
This is the reason why even 6 observations (6 samples) per site can be
too little!

### Normality

Often people say that the data should be normally distributed, if you
want to apply linear models.  
This is not exactly true: The normality should be met at each X value.
This is often difficult to check if the data does not have multiple
observations for each X and we don’t. In other words, the Y-data
contains the effects of all the explanatory variables (all metadata).
Therefore it is misleading to assess normality according the Y data.
Better choice is to apply a model and plot the pooled residuals. Since
the linear models should be LINEAR, the residuals should have the same
distance to fitted values (predicted means). If they do not, then the
model is not linear and the plot of residuals against fitted values will
show a curve or some other pattern. This means that the RESIDUALS should
be normally distributed, if linear models are used.

We know this about normal distribution:  
“It is a continuous probability distribution in which the random
variable can take any value.”  
With TSS-normalization, our values are “relative abundances”, so
something between 0 and 1.  
So relative abundances don’t meet the normality assumption and we cannot
use a model that assumes normality with TSS-normalized data.

Relative abundances are proportions, so the probability distribution to
use with them would be Beta distribution, but according to my
experience, Gamma distribution works also.  
Gamma distribution is used to model positive continuous variables that
represent time intervals between events, the Beta distribution is used
to model continuous variables that represent proportions or
probabilities.  
Gamma models work very well, I have had trouble with Beta

But what should be normally distributed even though it would be
calculated from relative abundances:  
Shannon’s diversity index.  
Some metadata variables (that are not %) would be normally distributed,
for instance pH.

### Fixed X

This assumption means that you know the values of X for each sample in
advance. You can have e.g. different fertilization rates, treatments,
etc. This assumption is violated, if you for instance study some
parameter of animals, which happen to be found somewhere and your X is
for example the age of animals. In other words, this assumption is
violated if the X is randomly selected.

We don’t have this problem.

### Zero trouble in Y-data

We probable have many zeros, at least in the ASV and genera levels

``` r
range(ASV_tableRA.mat)

plot(c(ASV_tableRA.mat))
```

And it seems we do. So if we would like to apply linear models using
taxa as the response variables, we might need to consider zero-inflated
models.

### Collinearity in X and relationships Y & X

If X-observations correlate, for instance weight and length, or water
depth and distance to the shoreline. If collinearity is ignored, one is
likely to end up with a confusing statistical analysis in which nothing
is significant, but where dropping one covariate can make the others
significant, or even change the sign of estimated parameters.  
collinearity can be detected with pairwise scatter plots comparing
covariates,

``` r
plot(MMB117metadata[, c(3:10)])
```

pH water and pH CaCl2 seem to correlate, as they should. But because of
the collinearity issue both of them cannot be used in a model.  
relationships: Y-values should be plotted against X-values to see which
X-values should be included in the model.

``` r
plot(MMB117metadata[, c(3:10, 20)])
```

This was a lazy persons solution. We can see that these is a
relationship with the diversity and the X-variables.  
And this also makes bilogically sense.

So if we would do linear models and use diversity as the response
variable, we might need to consider water %, pH (one of them) and SOM %

### Interactions?

This means that if we would do linear models, some of the X-variables
might have an interaction and this can be taken into account.  
Let’s try to visualize if we have something like this.  
I will categorize water-content

First summarize the moisture.

``` r
summary(MMB117metadata$Moisture)
```

Then we can categorize the moisture content and add them to the
metadata.

``` r
Water_FirstQ <- 38.50
Water_Median <- 42.16
Water_ThirdQ <- 54.97

MMB117metadata$MoistureCAT <- cut(MMB117metadata$Moisture, c(0, 45, 84.68))
```

Then we can check what we have.

``` r
MMB117metadata$MoistureCAT
levels(MMB117metadata$MoistureCAT)
```

And change the levels to something more informative.

``` r
levels(MMB117metadata$MoistureCAT) <- c("0 - 45 %", "45 - 84.7 %")
MMB117metadata$MoistureCAT
```

And then plot the pH and diversity

``` r
ggplot(subset(MMB117metadata, Site != "Control"), aes(x = pH_Ca, y = ASV_divSha, colour = Site)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    scale_color_manual(values = c("springgreen4", "chocolate4", "purple4", "royalblue")) +
    facet_grid(Site ~ MoistureCAT) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(colour = "black", fill = "white")) +
    theme(panel.border = element_rect(colour = "black")) +
    theme(legend.position = "none")
```

Looks like in our case we don’t have enough observations for this
interaction

## OPTIONAL: Possible contamination from control samples

We can have a closer look at the possible contaminants in our negative
control. We first identify all ASV that are present in our negative
control.  
Then we can check if these are present in our samples.

Extract only control samples.

``` r
ctrl_physeq <- physeq %>% subset_samples(Site=="Control")
```

Filter out ASVs that are present in the control (nmore than 1).

``` r
ctrl_taxa <- prune_taxa(taxa_sums(ctrl_physeq)>0, ctrl_physeq)
```

Get names of ASVs that are present in the control.

``` r
ctrl_taxa_names <- ctrl_taxa %>% taxa_names()
```

Extract only those taxa from the whole dataset.

``` r
ctrl_ASV_table <- prune_taxa(taxa_names(physeq) %in% ctrl_taxa_names,  physeq)
```

Plot a heatmap of the ASVs that are present in the control.

``` r
library(pheatmap)
pheatmap(otu_table(ctrl_ASV_table), cluster_rows = FALSE)
```

Square-root transformation of the data maybe looks better.

``` r
pheatmap(sqrt(otu_table(ctrl_ASV_table)), cluster_rows = FALSE)
```
