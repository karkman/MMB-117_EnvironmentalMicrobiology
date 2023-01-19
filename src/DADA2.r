.libPaths(c("~/projappl/project_rpackages_r421", .libPaths()))
libpath <- .libPaths()[1]

library(dada2)
library(phyloseq)
library(microViz)

# Paired-end Illumina data
 path <- "Work/Projects/naava/trimmed_data"
#list.files(path)

fnFs <- sort(list.files(path, pattern="_1_trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2_trimmed.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R2_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250, 200),
              maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(filtRs, err=errF, multithread=TRUE, pool=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

track <- cbind(out, rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample.names
track

#taxa50 <- assignTaxonomy(seqtab.nochim, "/Users/karkman/Work/Bioinfo/Databases/training_set.138_SSURef_NR99.fa.gz", minBoot=50, multithread=TRUE)
#taxa50 <- addSpecies(taxa50, "/Users/karkman/Work/Bioinfo/Databases/species_assignment.138_SSURef_NR99.fa.gz", tryRC = TRUE)

taxa80 <- assignTaxonomy(seqtab.nochim, "/Users/karkman/Work/Bioinfo/Databases/SILVA138/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
                        minBoot=80, multithread=TRUE)

# for pretty printing of the taxonomy
taxa.print <- taxa80
rownames(taxa.print) <- NULL
head(taxa.print)

naava <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE),
                 tax_table(taxa80))

dna <- Biostrings::DNAStringSet(taxa_names(naava))
names(dna) <- taxa_names(naava)
naava <- merge_phyloseq(naava, dna)

taxa_names(naava) <- paste0("ASV", seq(ntaxa(naava)))

## remove unwanted (both missing from data)
# mitochondria
naava <- prune_taxa(!(taxa_names(naava) %in% taxa_names(subset_taxa(naava, Family=="Mitochondria"))), naava)
# chloroplasts
naava <- prune_taxa(!(taxa_names(naava) %in% taxa_names(subset_taxa(naava, Order=="Chloroplast"))), naava)

naava_meta <- read.table("Work/Projects/naava/SraRunTable.txt", sep=",", header=TRUE, row.names=1)
naava_meta <- separate(naava_meta, Isolation_source, sep="_", into="SampleType", remove=FALSE)
sample_data(naava) <- sample_data(naava_meta)

saveRDS(naava, "Work/Projects/naava/naava.rds")
naava <- readRDS("Work/Projects/naava/naava.rds")

naava %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  comp_barplot(tax_level = "Genus", sample_order="default", n_taxa=19, facet_by = "SampleType") +
  coord_flip() +
  theme(axis.ticks.y=element_blank(),  axis.text.y=element_blank())

naava %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  comp_barplot(tax_level = "Genus", sample_order="default", n_taxa=19, facet_by = "SampleType") +
  coord_flip() +
  theme(axis.ticks.y=element_blank(),  axis.text.y=element_blank()) + facet_grid(Sample~.) +  facet_grid(
     rows = vars(DiseaseState),
     scales = "free", space = "free" # these options are critically important!
   )
