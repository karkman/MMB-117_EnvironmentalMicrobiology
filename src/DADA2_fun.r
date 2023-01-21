.libPaths(c("~/projappl/project_rpackages_r421", .libPaths()))
libpath <- .libPaths()[1]

library(dada2)
library(phyloseq)
library(microViz)

# Paired-end Illumina data
path <- "trimmed_data/fun_data"

fnFs <- sort(list.files(path, pattern="_1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2_trimmed.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R2_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, minLen=75,
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

taxa <- assignTaxonomy(seqtab.nochim, "/scratch/project_2007145/databases/sh_general_release_dynamic_29.11.2022.fasta",
                        minBoot=80, multithread=TRUE, tryRC = TRUE)

# for pretty printing of the taxonomy
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

physeq <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE),
                 tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(physeq))
names(dna) <- taxa_names(physeq)
physeq <- merge_phyloseq(physeq, dna)

taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))

physeq_meta <- read.table("doc/fun_meta.csv", sep=",", header=TRUE, row.names=1)
sample_data(physeq) <- sample_data(physeq_meta)

saveRDS(physeq, "outputs/physeq_fun.rds")
physeq <- readRDS("outputs/physeq_fun.rds")

physeq %>% ord_explore

physeq %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  comp_barplot(tax_level = "Genus", sample_order="default", n_taxa=19, facet_by = "SampleType") +
  coord_flip() +
  theme(axis.ticks.y=element_blank(),  axis.text.y=element_blank())

physeq %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  comp_barplot(tax_level = "Genus", sample_order="default", n_taxa=19, facet_by = "SampleType") +
  coord_flip() +
  theme(axis.ticks.y=element_blank(),  axis.text.y=element_blank()) + facet_grid(Sample~.) +  facet_grid(
     rows = vars(DiseaseState),
     scales = "free", space = "free" # these options are critically important!
   )
