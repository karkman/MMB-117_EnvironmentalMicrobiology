MMB-117
================

# Statistical analyses

During this excercise session we try to visualise our data and also test
our research questions.  
We will use the `phyloseq` and `microViz` packages to visualise the data
and `vegan` package to test the hypothesis.  
Below you find some examples, but you need to modify the code to
visualise and test all possible combinations (that make sense).

## Setup

``` r
.libPaths(c("/projappl/project_2013123/project_rpackages_r421", .libPaths()))
libpath <- .libPaths()[1]
library(tidyverse)
library(phyloseq)
library(vegan)
```

First we need to read in the data and make the the objects for plotting.

Read in the phyloseq object and extract the metadata from the object.  
Also add the Shannon diversity index to the metadata.

``` r
physeq <- readRDS("physeq.rds")
physeq <- subset_samples(physeq, Site != "Control")

MMB117metadata <- sample_data(physeq) %>%
    as_tibble() %>%
    mutate(Sample = rownames(sample_data(physeq))) %>%
    as.data.frame()

MMB117metadata$ASV_divSha <- diversity(otu_table(physeq), index = "shannon")

View(MMB117metadata)
```

## Transformations/normalisations

We can use phyloseq to normalise and rarefy the data. And use a fucntion
from vegan to do the CLR transformation.

Relative abundances.

``` r
physeq_ra <- transform_sample_counts(physeq, function(x) x / sum(x))
```

For CLR transformation we use `decostand` function from vegan package.  
We need to extract the ASV table from the phyloseq object and then make
a new objet and store the transformed count to that object.  
But you can find many other packages that can do this.

``` r
ASV_table <- physeq %>% otu_table()
ASV_table <- decostand(ASV_table, "clr", 1, pseudocount=1)
physeq_clr <- physeq
otu_table(physeq_clr) <- otu_table(ASV_table, taxa_are_rows=FALSE)
```

Rarefaction of the data. We can use the `rarefy_even_depth` function
from phyloseq package.

``` r
physeq_rare <- rarefy_even_depth(physeq, sample.size=min(sample_sums(physeq)), rngseed=1)
```

Extract the differntly normalised data for further analyses.  
Raw counts as an example.

``` r
ASV_table <- physeq %>%
    otu_table() %>%
    as.data.frame()
```

Outcomes of the data exploration: - No reason to remove any values
(outliers) - Heteroscedasticity in all variables except pH_Ca and
Shannon diversity index this is within acceptable levels - Normality: If
we use relative abundances, they are not normally distributed. Shannon’s
diversity index would be - Fixed X: we don’t have this problem -
Collinearity in X and Relationships Y & X: We need to select only one pH
and looks like there are relationships between Y-data and all X-values.
So perhaps we should model the influence of pH_Ca, SOM and GWC on
communities. - Interactions: we don’t have enough observations for a
model with interactions. - About GWC: we can include it in the analysis
for practice. But actually there is an issue and it is the amount of
snow. We didn’t have a system for removing the snow or taking the amount
of snow into account. So GWC might not be that interesting although it
could explain a lot “numerically”.

Biostatistical analyses: One of the original research questions was to
see how the diversity & community structure is influenced by human
activities (different sites) hypothesis was that the diversity would be
lower in the gas station.  
Based on data exploration, we also should ask what is the influence of
pH_Ca, SOM and GWC on communities. We can answer these questions with
ordination analyses.  
We can model pH_Ca and diversity as linear vectors and SOM % and GWC %
as nonlinear surfaces (GAM) test that NMDS works and check the
dimensions.

Now we have rarefied data, so we could count the species richness.
THere’s a function `specnumber` in vegan package that can be used for
this.  
Below you have an example how to use it. How would you add the species
richness to the metadata?

``` r
specnumber(otu_table(physeq_rare))
```

We can make a barplot of the most abundant ASVs in our data using
functions from `microViz` package.  
Read the [microViz
documentation](https://david-barnett.github.io/microViz/) to see how to
do this.

``` r
library(microViz)

# barplot code here
```

Then we can make a PCoA plot of the data. We can use the `ordinate` and
`plot_ordination` functions from phyloseq package.  
Make different plots for different normalisations/transformations. And
change the method and distance accordingly.  
The plot functions uses ggplot2 package, so you can modify the plot as
you like.

Example with raw counts.

``` r
pcoa <- ordinate(physeq, method="PCoA", distance="bray")
plot_ordination(physeq, pcoa, color="Site")
```

Then we can make a bit more sophisticated NMDS plots.  
First we need to make the NMDS with `metaMDS` function from vegan
package.

``` r
plot(metaMDS(ASV_table, distance="bray", k=2), type="text", display="sites")
```

Save the NMDS

``` r
MMB117_NMDS<-metaMDS(ASV_table, distance="bray", k=2)
```

Add colors to metadata

``` r
levels(MMB117metadata$Site)
MMB117metadata$color<-rep(1, nrow(MMB117metadata))
MMB117metadata<- within(MMB117metadata, color[Site=="Field"]<-"greenyellow")
MMB117metadata<- within(MMB117metadata, color[Site=="Forest"]<-"forestgreen")
MMB117metadata<- within(MMB117metadata, color[Site=="Gas_station"]<-"darkslategray4")
MMB117metadata<- within(MMB117metadata, color[Site=="Park"]<-"darkkhaki")
```

Environmental fitting of the NMDS with diversity and pH

``` r
MMB117EF<-envfit(MMB117_NMDS ~ ASV_divSha + pH_Ca, MMB117metadata, permutations=999)
MMB117EF
```

Plot the NMDS with diversity and pH

``` r
NMDSoplot<-ordiplot(MMB117_NMDS,type="n", xlim=c(-1.8,1.8),
                  ylim=c(-1.8,1.8),cex.axis = 1.5, cex.lab = 1.5)
with(MMB117metadata, points(MMB117_NMDS$points,pch=15, cex=2, col=MMB117metadata$color))
plot(MMB117EF, col="cyan4", cex=1)
with(MMB117metadata, ordisurf(MMB117_NMDS,SOM, add = TRUE, col = "grey20", cex=2))
#with(MMB117metadata_noneg, ordisurf(MMB117_NMDS,Moisture, add = TRUE, col = "royalblue"))
identify(NMDSoplot, "sites", labels = MMB117metadata$Sample, cex=1)
with(MMB117metadata, legend(1.2,1.8, legend= levels(Site), cex=1.5, bty= "n",
       col=c("greenyellow", "forestgreen", "darkslategray4","darkkhaki"), pch=c(15,15,15,15)))
legend(1.25,0.98,"SOM %",cex=1.5,lty=1,col="grey20",bty= "n")
```

Identify samples by clicking them. I would click only those that are
“outliers”, to get the idea how SOM influences on communities Hit “esc”
when you are ready. Check that your cursor is in the terminal if nothing
happens.

Statistical interpretation:

``` r
NMDS_SOM_surf <- ordisurf(MMB117_NMDS ~ SOM, MMB117metadata)
summary(NMDS_SOM_surf)
```

ordisurf fits a GAM model and accepts nonlinear variables We have an
intercept model and the estimate is a mean of our response variable (SOM
%) R-sq.(adj) = 0.788 suggests that ~79 % of the variance is explained
by SOM % (pretty good!) Deviance explained = 86.4% indicates the
goodness of fit, which is also very good in our case. So what the model
suggests is that SOM % explains ~79 % of the variance in the community
structure (=beta diversity) (p \< 0.05 ) Interestingly, the diversity is
growing to the direction of the gas station! Why? It is lowest in the
forest and field sites Why?? 6:41

Same plot with GWC including pH_Ca and div as linear vectors and surface
fitting GWC %

``` r
NMDSoplot<-ordiplot(MMB117_NMDS,type="n", xlim=c(-1.8,1.8),
                    ylim=c(-1.8,1.8),cex.axis = 1.5, cex.lab = 1.5)
with(MMB117metadata, points(MMB117_NMDS$points,pch=15, cex=2, col=MMB117metadata$color))
plot(MMB117EF, col="cyan4", cex=1)
#with(MMB117metadata_noneg, ordisurf(MMB117_NMDS,SOM, add = TRUE, col = "grey20", cex=2))
with(MMB117metadata, ordisurf(MMB117_NMDS,Moisture, add = TRUE, col = "royalblue", cex=2))
identify(NMDSoplot, "sites", labels = MMB117metadata$Sample, cex=1)
with(MMB117metadata, legend(1.2,1.8, legend= levels(Site), cex=1.5, bty= "n",
                                  col=c("greenyellow", "forestgreen", "darkslategray4","darkkhaki"), pch=c(15,15,15,15)))
legend(1.25,0.98,"GWC %",cex=1.5,lty=1,col="royalblue",bty= "n")
```

Identify samples by clicking them. I would click only those that are
“outliers”, to get the idea how GWC influences on communities.  
Hit “esc” when you are ready. Check that your cursor is in the terminal
if nothing happens.

Do the statistical interpretation for GWC yourself.

We can also model influence of the variables on metadata with Permanova:

``` r
adonis2(ASV_table ~ Site + pH_Ca + SOM + Moisture, data=MMB117metadata, permutations=9999, by = "terms", na.action = na.omit)

adonis2(ASV_table ~ pH_Ca + SOM + Moisture, data=MMB117metadata, permutations=9999, by = "terms", na.action = na.omit)
```

Now SOM is not significant! Let’s see what happens if we leave away GWC
(moisture).

``` r
adonis2(ASV_table ~ pH_Ca + SOM, data=MMB117metadata, permutations=9999, by = "terms", na.action = na.omit)
```

So according to permanova SOM is not significant.  
What happens if we leave out pH?

``` r
adonis2(ASV_table ~  SOM, data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)
```

What if we include pH and GWC?

``` r
adonis2(ASV_table ~ pH_Ca + Moisture, data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)
```

So maybe we should include only site in our permanova, or then pH
without the site.  
Why is this and we needed to check how the results change?

``` r
adonis2(ASV_table ~ pH_Ca , data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)
adonis2(ASV_table ~ Site, data=MMB117metadata_noneg, permutations=9999, by = "terms",na.action = na.omit)
```

Which one we should have, Site or pH?

Influence of the site on alpha diversity:  
We can compare the diversities of the sites with for instance a t-test
(remember the normality):

``` r
library(ggpubr)
ggplot(MMB117metadata, aes(x=Site, y=ASV_divSha)) +
  geom_point(aes(fill=factor(Site)), size=3, shape=21, colour="grey20",alpha=0.7,
             position=position_jitter(width=0.01, height=0.01)) +
  geom_boxplot (outlier.colour = NA, fill=NA, colour="grey20") +
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(colour = "black", fill = "white"))+
  theme(axis.text.x = element_text(angle = 35,hjust=0.7,vjust=0.8))+
  theme(panel.border = element_rect(colour = "black"))+
  #theme(axis.text.x=element_blank())+
  scale_fill_manual(values=c("greenyellow", "forestgreen", "darkslategray4","darkkhaki"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  theme(legend.position="none")+
  geom_hline(yintercept = mean(na.omit(MMB117metadata$ASV_divSha)), linetype = 2)+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "Forest", hide.ns = TRUE)
```

So significantly higher in Gas station and park than in Forest.  
What does “significantly” mean? We can also get the comparison in a
table format like this:

``` r
compare_means(ASV_divSha ~ Site,  data = MMB117metadata,
              ref.group = "Forest", method = "t.test")
```

**So what you will write in your reports?**
