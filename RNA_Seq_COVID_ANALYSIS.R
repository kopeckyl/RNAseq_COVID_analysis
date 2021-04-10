### RNA-seq analysis for COVID and SARS expression comparison

## Analysis and code by Lucas K. Bobadilla
# If use for publication please cite: 


# Material and methods


## Data source and pre-processing

#Data for this comparison were extracted from two GEO datasets (GSE147507 and GSE122876).
#Both datasets were single read RNA-seq files from Calu-3 cells infected with the virus after 24h.
#Both experiments were conducted using triplicates for each condition.
#Raw fastq files from both datasets were downloaded from NCBI using the SRA Toolkit and
#data quality was accessed using fastQC reports.
#All raw fastq files were trimmed using
#Trimmomatic set for LEADING:28, TRAILING:28, SLIDINGWINDOW:4:30 and MINLEN:50.
# After trimming low quality reads, libraries we pseudo-aligned to the homo sapiens genome GRCh38 using
# Kallisto without bootstrapping to get the transcript abundance.



# Analysis

## Packages

library(sva)
library(edgeR)
library(tximport)
library(tidyverse)
library(biomaRt)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggVennDiagram)
library(goseq)
library(kableExtra)




## metadata


# sample information - meta data
samples <-read.table("samples.txt", header=TRUE)
samples$condition <- factor(samples$condition, levels = c("Control","COVID19", "MERS"))
samples$batch <- factor(samples$batch)
rownames(samples) <- samples$sample
samples <- samples[,-1]

kable(samples)


## Load biomart query


#biomart query

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
tx2gene <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id","hgnc_symbol", "chromosome_name"), mart = ensembl)




## Load Abundance reads

#get abundance files path
files <- file.path("data/all_data/", rownames(samples), "abundance.tsv")
names(files) <- rownames(samples)
files

# Load files using tximport - transcript level
txi.kallisto.covid <- tximport(files, type = "kallisto", tx2gene = tx2gene, txOut = TRUE)

#change to gene level

gene.kallisto.covid <- summarizeToGene(txi.kallisto.covid,tx2gene)



## Apply model matrix and Combat_seq function

# matrix design
full.mat <- model.matrix(~condition, data = samples)

null.mat <- model.matrix(~1, data = samples)

#solving batch effect
batch <- samples$batch
group <- samples$condition
adjusted <- ComBat_seq(gene.kallisto.covid$counts, batch=batch, group=group)

# transforming data to integer
adj <- sapply(data.frame(adjusted), as.integer)
rownames(adj) <- rownames(adjusted)


dds <- DESeqDataSetFromMatrix(countData = adj,
                              colData = samples,
                              design = ~ condition)

rld <- rlog(dds)
plotPCA(rld, intgroup = "condition")
#ggsave(file = "PCA.jpeg", dpi = 300, units = "cm", height = 10, width = 24)


## Apply EdgeR model and plot dispersion

dgList <- DGEList(counts=adjusted, genes=rownames(adjusted), group=samples$condition)

countsPerMillion <- cpm(dgList) #get counts
countCheck <- countsPerMillion > 1 # create a conditional matrix for values > 1
keep <- which(rowSums(countCheck) >= 1) #filter values
dgList <- dgList[keep,]

dgList_norm <- calcNormFactors(dgList, method="TMM")

#dispersion
disp <- estimateDisp(dgList_norm, design=full.mat, robust=TRUE)

## plot dispersions
plotBCV(disp)





## volcano plot and DEG for COVID

# DEG analysis

fit <- glmFit(disp, full.mat)
lrt <- glmLRT(fit, coef=2) # fit only for covid vs control

DEG_genes_covid <- topTags(lrt, adjust.method = "BH", p.value = 0.001, n = 400) # get genes significant at FDR < 0.001


# get all genes results
result_EDGER_covid <- topTags(lrt, n=nrow(lrt$table), adjust.method = "BH")

# get gene IDs from transcripts
FDR_genes_covid <- result_EDGER_covid$table # convert to dataframe

FDR_genes_significant_covid <- DEG_genes_covid$table # convert to dataframe

plotMD(lrt) # all significant ones

FDR_genes_covid %>%
  ggplot(aes(y = -log10(FDR), x = logFC)) +
  geom_point() +
  geom_point(aes(y = -log10(FDR), x = logFC), col ="red", data = FDR_genes_significant_covid) + theme_light() +
  labs(title = "Volcano plot - COVID vs Control")
#ggsave("volcano_COVID.jpeg", dpi = 300)





## Volcano plot for MERS

# DEG analysis


lrt <- glmLRT(fit, coef=3) # fit model for MERS vs control
DEG_genes_mers <- topTags(lrt, adjust.method = "BH", p.value = 0.001, n = 400)
result_EDGER_mers <- topTags(lrt, n=nrow(lrt$table), adjust.method = "BH")

# get gene IDs from transcripts
FDR_genes_mers <- result_EDGER_mers$table

FDR_genes_significant_mers <- DEG_genes_mers$table

#plotMD(lrt) # all significant ones

FDR_genes_mers %>%
  ggplot(aes(y = -log10(FDR), x = logFC)) +
  geom_point() +
  geom_point(aes(y = -log10(FDR), x = logFC), col ="red", data = FDR_genes_significant_mers) + theme_light() +
  labs(title = "Volcano plot - MERS vs Control")

#ggsave("volcano_MERS.jpeg", dpi = 300)




# Find genes in common between comparison with control

# find genes in top common genes between disease when compared to control
common_disease <- FDR_genes_significant_covid %>%
  dplyr::select(genes, logFC_covid = logFC,FDR_covid = FDR) %>%
  inner_join(FDR_genes_significant_mers %>%
               dplyr::select(genes, logFC_mers= logFC,FDR_mers = FDR),
             by ="genes")

write.csv(common_disease, "top100_common_genes.csv")

common_disease_tidy <- common_disease %>%
  dplyr::select(genes, covid=FDR_covid,mers=FDR_mers) %>%
  gather(covid:mers, key = "condition", value = "FDR") %>%
  left_join(common_disease %>%
              dplyr::select(genes, covid=logFC_covid,mers=logFC_mers) %>%
              gather(covid:mers, key = "condition", value = "logFC"))


common_genes <- common_disease %>%
  dplyr::select(ensembl_gene_id=genes,everything() ) %>%
  inner_join(tx2gene %>% dplyr::select(-ensembl_transcript_id_version)) %>%
  distinct() %>%
  dplyr::select(hgnc_symbol, ensembl_gene_id,
                chromosome_name, logFC_covid,
                logFC_mers) %>%
  arrange(chromosome_name, desc(logFC_covid))


all_genes <- FDR_genes_significant_covid %>%
  mutate(condition = "COVID") %>%
  full_join(FDR_genes_significant_covid %>%
              mutate(condition = "MERS"))

#write_csv(all_genes,"top400_genes.csv")

# get all significant genes that were common

common_disease_all <- FDR_genes_covid %>%
  filter(FDR <=0.001) %>%
  dplyr::select(genes, logFC_covid = logFC,FDR_covid = FDR) %>%
  inner_join(FDR_genes_mers %>%
               filter(FDR <=0.001) %>%
               dplyr::select(genes, logFC_mers= logFC,FDR_mers = FDR),
             by ="genes")


common_disease_all <- common_disease_all %>%
  dplyr::select(ensembl_gene_id=genes,everything() ) %>%
  inner_join(tx2gene %>%
               dplyr::select(-ensembl_transcript_id_version)) %>%
  distinct() %>%
  dplyr::select(hgnc_symbol, ensembl_gene_id,
                chromosome_name, logFC_covid,
                logFC_mers) %>%
  arrange(chromosome_name, desc(logFC_covid))

#write_csv(common_disease_all, "common_genes.csv")


#cytokine genes plot
common_disease_all %>% arrange(desc(logFC_covid), desc(logFC_mers)) %>%
  filter(str_detect(hgnc_symbol, "^IL")) %>%
  gather(logFC_covid:logFC_mers, key = "Disease", value = "LogFC") %>%
  separate(Disease, into = c("disc", "Disease"), sep = "_") %>%
  ggplot(aes(x = hgnc_symbol, y = LogFC, fill = Disease)) +
  geom_col(position = position_dodge(.5), col = "black",width = 0.5) +
  geom_hline(yintercept = 0) +
  theme_light() +
  scale_fill_brewer(palette="Dark2") +
  labs(x = "Gene", title = "Differentially expressed genes involved with cytokine response")

#ggsave("IL_genes.jpeg", dpi = 300, units = "cm", width = 24, height = 12)

#GBP genes
common_disease_all %>% arrange(desc(logFC_covid), desc(logFC_mers)) %>%
  filter(str_detect(hgnc_symbol, "^GBP")) %>%
  gather(logFC_covid:logFC_mers, key = "Disease", value = "LogFC") %>%
  separate(Disease, into = c("disc", "Disease"), sep = "_") %>%
  ggplot(aes(x = hgnc_symbol, y = LogFC, fill = Disease)) +
  geom_col(position = position_dodge(.5), col = "black",width = 0.5) +
  geom_hline(yintercept = 0) +
  theme_light() +
  scale_fill_brewer(palette="Dark2") +
  labs(x = "Gene", title = "Differentially expressed Guanylate-binding proteins (GBP) genes")

#ggsave("GBP_genes.jpeg", dpi = 300, units = "cm", width = 24, height = 12)


#TNF genes

common_disease_all %>% arrange(desc(logFC_covid), desc(logFC_mers)) %>%
  filter(str_detect(hgnc_symbol, "^TNF")) %>%
  gather(logFC_covid:logFC_mers, key = "Disease", value = "LogFC") %>%
  separate(Disease, into = c("disc", "Disease"), sep = "_") %>%
  ggplot(aes(x = hgnc_symbol, y = LogFC, fill = Disease)) +
  geom_col(position = position_dodge(.5), col = "black",width = 0.5) +
  geom_hline(yintercept = 0) +
  theme_light() +
  scale_fill_brewer(palette="Dark2") +
  labs(x = "Gene", title = "Differentially expressed tumour necrosis factor (TNF) genes")

#ggsave("TNF_genes.jpeg", dpi = 300, units = "cm", width = 24, height = 12)



## Venn diagram
venn_mers <- FDR_genes_mers %>%
  filter(FDR <=0.001)

rownames(venn_mers) <- venn_mers$genes

venn_covid <- FDR_genes_covid %>%
  filter(FDR <=0.001)

rownames(venn_covid) <- venn_covid$genes

venn_plot <- list(COVID = rownames(venn_covid),
                  MERS = rownames(venn_mers))


ggVennDiagram(venn_plot)
#ggsave(file = "venn_diagram.jpeg", dpi = 300, height = 10, width = 20, units = "cm")




## list of genes up and down regulated

# covid
FDR_genes_significant_covid %>%
  mutate(DEG = if_else(logFC < 0, "Down-regulated", "Up-regulated")) %>%
  group_by(DEG) %>%
  summarize(DEG_count = n())

FDR_genes_significant_mers %>%
  mutate(DEG = if_else(logFC < 0, "Down-regulated", "Up-regulated")) %>%
  group_by(DEG) %>%
  summarize(DEG_count = n())




## GO term annotation

# get GO terms
query <- common_disease_all %>%
  mutate(DEG_covid = if_else(logFC_covid < 0, "Down-regulated", "Up-regulated"),
         DEG_mers = if_else(logFC_covid < 0, "Down-regulated", "Up-regulated")) %>%
  dplyr::select(hgnc_symbol, DEG_covid, DEG_mers,chromosome_name)

common_disease_all %>%
  arrange(desc(logFC_covid),desc(logFC_mers))


query_result <- getBM(attributes = c("hgnc_symbol", "go_id",
                                     "name_1006","definition_1006","namespace_1003"),
                      filters = "hgnc_symbol",
                      values = query$hgnc_symbol,
                      mart = ensembl)


Enrichment_analysis <- query %>%
  left_join(query_result) %>%
  dplyr::select(-chromosome_name) %>%
  filter(namespace_1003 == "biological_process") %>%
  dplyr::select(-namespace_1003)

colnames(Enrichment_analysis) <- c("gene", "COVID", "MERS",
                                   "go_ID", "Biological_process",
                                   "Biological_description")

write.csv(Enrichment_analysis,"Enrichment_analysis.csv")


## GO-seq enrichment analysis
library("org.Hs.eg.db")
deg_genes <- as.vector(common_disease_all$ensembl_gene_id)
measure_genes <- as.vector(FDR_genes_mers$genes)

gene_vector <- as.integer(measure_genes%in%deg_genes)
names(gene_vector) <- measure_genes


pwf <- nullp(gene_vector,"hg19","ensGene")
head(pwf)

GO.wall <- goseq(pwf,"hg19","ensGene",test.cats=c("GO:BP"))

GO.wall %>%
  top_n(30, wt=-over_represented_pvalue) %>%
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc,
             y=term,
             colour=over_represented_pvalue,
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p-value", size="Count", title = "GO terms - top 30")

#ggsave("GO-terms.jpeg", dpi = 300)



