# Nature Scientific Data Paper
# Author: Kelsey Montgomery
# Date: 20190604

# Directory for figure output
dir.create("files")

# Process functions
downloadFile_version <- function(id , version){
  fread(synGet(id, version = version)$path, data.table = F)
}

# Dependencies
library(spatstat)
library(ggrepel)
library(cowplot)
library(synapser)
library(dplyr)
library(biomaRt)
library(data.table)
library(ggplot2)
library(plyr)
library(CovariateAnalysis)
library(edgeR)
library(stringr)

# Login to Synapse - https://r-docs.synapse.org/articles/manageSynapseCredentials.html
synLogin()

# Download counts (DLPFC - MSSM)
COUNT_ID = 'syn17346208'
count = downloadFile_version(COUNT_ID, version = 2) %>% data.frame()
count$transcript_id.s. = NULL

# Download gene lengths (DLPFC - MSSM)
genelen_CMC = downloadFile_version('syn17346397', version = 2) %>%
  tidyr::gather(sampleID, Length, -gene_id, -`transcript_id(s)`) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(Length = median(Length, na.rm = T)) %>%
  ungroup() %>% data.frame()

# Download counts (DLPFC - HBCC)
COUNT_ID = 'syn17894685'
count_HBCC = downloadFile_version(COUNT_ID, version = 4) %>% data.frame()
count_HBCC$transcript_id.s. = NULL

# Download gene lengths (DLPFC)
genelen_HBCC = downloadFile_version('syn18324060', version = 3) %>% 
  tidyr::gather(sampleID, Length, -gene_id, -`transcript_id.s.`) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(Length = median(Length, na.rm = T)) %>%
  ungroup() %>% data.frame()

# Download clinical metadata 
CLINICAL_ID = 'syn3354385'
clinical = downloadFile_version(CLINICAL_ID, version = 4)

# Download RNASeq metadata
METADATA_QC_DLPFC_ID = 'syn18358379' 
metadata = downloadFile_version(METADATA_QC_DLPFC_ID, version = 3)

# Join clinical and RNASeq metadata 
md = right_join(clinical, metadata, by = c("Individual ID" = "Individual_ID")) %>% 
  dplyr::mutate(Dx = forcats::fct_recode(Dx, AFF_BP = "BP", AFF_BP = "AFF", Other = "undetermined", Control = "Control", SCZ = "SCZ"))

# Compute read pair metrics and add Institution-Dx variable
md <- md %>%
  mutate(MappedRead_Pairs = Mapped_Reads/2) %>%
  mutate(`Institution-Dx` = paste0(`Institution`, "-", `Dx`)) %>% 
  mutate(TotalRead_Pairs = Total_Reads/2)

# Join HBCC and MSSM counts
count = full_join(count, count_HBCC, by = c("gene_id"))

# Bind gene lengths
genelen = bind_rows(genelen_CMC, genelen_HBCC) %>%
  dplyr::group_by(gene_id) %>% 
  dplyr::summarize(Median = median(Length)) %>% 
  dplyr::rename(Length = Median)

# Get GC content from biomart
backgroundGenes = data.frame(gene_id = count$gene_id) %>%
  dplyr::mutate(id = gene_id) %>%
  tidyr::separate(id, c('ensembl_gene_id','position'), sep = '\\.')

# Define biomart object
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "jan2019.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")

# Query biomart
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "percentage_gene_gc_content", "gene_biotype", "chromosome_name"),
                       filters = "ensembl_gene_id", values = backgroundGenes$ensembl_gene_id,
                       mart = mart)
# Get GC content
gc_content = Ensemble2HGNC %>%
  dplyr::left_join(backgroundGenes) %>% 
  dplyr::select(gene_id, percentage_gene_gc_content, chromosome_name) %>%
  unique
rownames(gc_content) = gc_content$gene_id

# Sex check with XIST and UTY expression - ENSG00000229807.10 and ENSG00000183878.15
sex_counts = gc_content %>% 
  left_join(count) %>%
  dplyr::select(-one_of("percentage_gene_gc_content")) %>%
  filter(chromosome_name == "X" |chromosome_name == "Y") %>% 
  tidyr::gather(key = item, value = value, -c(gene_id, chromosome_name)) %>%
  mutate(value = log(value)) %>%
  dplyr::rename(`counts(log)`= value) %>% 
  dplyr::rename(SampleID = item) %>%
  left_join(md[,c("SampleID", "Sex", "Institution")])

filt <- sex_counts %>% 
  filter(gene_id == "ENSG00000229807.10" | gene_id == "ENSG00000183878.15") %>% 
  dplyr::select(-one_of("chromosome_name")) %>% 
  tidyr::spread(key = gene_id, value = `counts(log)`) %>% 
  mutate(XIST = as.numeric(`ENSG00000229807.10`)) %>% 
  mutate(UTY = as.numeric(`ENSG00000183878.15`)) %>% 
  mutate(UTY = ifelse(UTY == -Inf, 0, UTY)) %>% 
  mutate(XIST = ifelse(XIST == -Inf, 0, XIST)) %>% 
  mutate(Institution = ifelse(Institution == "NIMH-HBCC", "NIMH_HBCC", Institution))


# Covariate relationships
FactorCovariates <- c('Individual ID', "Institution", "Reported Gender", "Library_Batch", "Ribozero_Batch", "Flowcell_Batch", "Dx", "Institution-Dx")
ContCovariates <- c("Age of Death", "PMI (in hours)", "RIN", "Mapped_Reads", "MappedRead_Pairs", "Intragenic_Rate", "Intronic_Rate", "Intergenic_Rate",
                    "Genes_Detected", "Expression_Profiling_Efficiency", "rRNA_Rate", "Total_Reads", "TotalRead_Pairs","Percent_Aligned")
# Ages over 90
md <- mutate(md, `Age of Death` = ifelse(`Age of Death` == "90+", "90", `Age of Death`))

# Find inter relation between factor covariates
covariates = md[,c(FactorCovariates,ContCovariates),drop=F]
covariates[,FactorCovariates] <- data.frame(lapply(covariates[,FactorCovariates],function(x){
  x <- sapply(x,function(y){str_replace_all(as.character(y),'[^[:alnum:]]','_')})}))
rownames(covariates) <- md$SampleID

# Convert factor covariates to factors
covariates[,FactorCovariates] = lapply(covariates[,FactorCovariates], factor)
covariates[,ContCovariates] = lapply(covariates[,ContCovariates], function(x){
  x = as.numeric(as.character(gsub('[\\,\\%]','',x)))
})

#Summary statistics 
stat <- as_tibble(list(Covariates = c("RIN", "Intergenic Rate", "Total Read Pairs", "Intronic Rate", "Mapped Read Pairs", "rRNA Rate"),
                  Mean = c(mean(covariates$RIN), mean(covariates$Intergenic_Rate), mean(covariates$TotalRead_Pairs), 
                           mean(covariates$Intronic_Rate), mean(covariates$MappedRead_Pairs), mean(covariates$rRNA_Rate)),
                  SD = c(sd(covariates$RIN), sd(covariates$Intergenic_Rate), sd(covariates$TotalRead_Pairs), 
                         sd(covariates$Intronic_Rate), sd(covariates$MappedRead_Pairs), sd(covariates$rRNA_Rate))))
stat

# Figure 1
my.theme <- theme_bw(12) %+replace% theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# RIN
p = list()
p[[1]] = ggplot(covariates, aes(x = Dx, y = RIN)) + geom_boxplot()
p[[1]] = p[[1]] + ggtitle('RIN') + my.theme + xlab("Diagnosis")

# Intergenic Rate
p[[2]] = ggplot(covariates, aes(x = Dx, y = Intergenic_Rate)) + geom_boxplot()
p[[2]] = p[[2]] + ggtitle('Intergenic Rate') + my.theme + ylab("Intergenic Rate") + xlab("Diagnosis")

# Intronic Rate
p[[3]] = ggplot(covariates, aes(x = Dx, y = Intronic_Rate)) + geom_boxplot()
p[[3]] = p[[3]] + ggtitle('Intronic Rate') + my.theme + ylab("Intronic Rate") + xlab("Diagnosis")

# Mapped Reads
p[[4]] = ggplot(covariates, aes(x = Dx, y = MappedRead_Pairs)) + geom_boxplot()
p[[4]] = p[[4]] + ggtitle('Mapped Read Pairs') + my.theme + ylab("Mapped Read Pairs") + xlab("Diagnosis")

# Total Reads
p[[5]] = ggplot(covariates, aes(x = Dx, y = TotalRead_Pairs)) + geom_boxplot()
p[[5]] = p[[5]] + ggtitle('Total Read Pairs') + my.theme + ylab("Total Read Pairs") + xlab("Diagnosis")

# rRNARate
p[[6]] = ggplot(covariates, aes(x = Dx, y = rRNA_Rate)) + geom_boxplot()
p[[6]] = p[[6]] + ggtitle('rRNA Rate') + my.theme + ylab("rRNA Rate") + xlab("Diagnosis")

plotVar = plot_grid(p[[1]],p[[2]],p[[3]], p[[4]], p[[5]], p[[6]], align = "h", ncol = 3, labels = c("a","b","c","d","e","f"))
save_plot("./files/figure1.pdf", plotVar,
          ncol = 3, 
          nrow = 2, 
          base_aspect_ratio = 1
)

# Supplementary Figure 1
my.theme <- theme_bw() %+replace% theme(legend.position = 'top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.title = element_text(hjust = 0.5))

# RIN
p = list()
p[[1]] = ggplot(covariates, aes(x = `Institution-Dx`, y = RIN)) + geom_boxplot()
p[[1]] = p[[1]] + ggtitle('RIN') + my.theme + xlab("Institution-Diagnosis")

# Intergenic Rate
p[[2]] = ggplot(covariates, aes(x = `Institution-Dx`, y = Intergenic_Rate)) + geom_boxplot()
p[[2]] = p[[2]] + ggtitle('Intergenic Rate') + my.theme + ylab("Intergenic Rate") + xlab("Institution-Diagnosis")

# Intronic Rate
p[[3]] = ggplot(covariates, aes(x = `Institution-Dx`, y = Intronic_Rate)) + geom_boxplot()
p[[3]] = p[[3]] + ggtitle('Intronic Rate') + my.theme + ylab("Intronic Rate") + xlab("Institution-Diagnosis")

# Mapped Reads
p[[4]] = ggplot(covariates, aes(x = `Institution-Dx`, y = MappedRead_Pairs)) + geom_boxplot()
p[[4]] = p[[4]] + ggtitle('Mapped Read Pairs') + my.theme + ylab("Mapped Read Pairs") + xlab("Institution-Diagnosis")

# Total Reads
p[[5]] = ggplot(covariates, aes(x = `Institution-Dx`, y = TotalRead_Pairs)) + geom_boxplot()
p[[5]] = p[[5]] + ggtitle('Total Read Pairs') + my.theme + ylab("Total Read Pairs") + xlab("Institution-Diagnosis")

# rRNARate
p[[6]] = ggplot(covariates, aes(x = `Institution-Dx`, y = rRNA_Rate)) + geom_boxplot()
p[[6]] = p[[6]] + ggtitle('rRNA Rate') + my.theme + ylab("rRNA Rate") + xlab("Institution-Diagnosis")

plotSupplementary = plot_grid(p[[1]],p[[2]],p[[3]], p[[4]], p[[5]], p[[6]], align = "h", ncol = 3, labels = c("a","b","c","d","e","f"))
save_plot("./files/supplementaryFigure1.pdf", plotSupplementary,
          ncol = 3, 
          nrow = 2, 
          base_aspect_ratio = 1
)


# Filter genes
# Remove genes that have less than 1 cpm counts in at least 50% of samples per Dx and per Dx.Reported Gender. Also remove genes with missing gene length and percentage GC content.
rownames(count) = count$gene_id
count$gene_id = NULL

genesToAnalyze = plyr::dlply(md, .(Dx, `Reported Gender`), .fun = function(mtd, count){
  processed.counts = getGeneFilteredGeneExprMatrix(count[,mtd$SampleID],
                                                   MIN_GENE_CPM=1, 
                                                   MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.5)
  processed.counts$filteredExprMatrix$genes
}, count)

genesToAnalyze = unlist(genesToAnalyze) %>% 
  unique() %>% 
  intersect(gc_content$gene_id[!is.na(gc_content$percentage_gene_gc_content)]) %>% 
  intersect(genelen$gene_id[!is.na(genelen$Length)])

processed_counts = getGeneFilteredGeneExprMatrix(count[genesToAnalyze, ], MIN_GENE_CPM=0, MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0)
pct.pc = processed_counts$filteredExprMatrix$genes %>%
  dplyr::rename(gene_id = genes) %>%
  left_join(backgroundGenes) %>%
  left_join(Ensemble2HGNC) %>%
  dplyr::group_by(gene_biotype) %>%
  dplyr::summarise(fraction  = n()) %>%
  dplyr::filter(fraction > 100) %>%
  dplyr::mutate(fraction = fraction/length(processed_counts$filteredExprMatrix$genes[,1]))

PC <- prcomp(voom(processed_counts$filteredExprMatrix$counts)$E, scale.=T, center = T)

# Plot first 2 PCs
plotdata <- data.frame(SampleID=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])

# Percentage from each PC
eigen <- PC$sdev^2
pc1 <- eigen[1]/sum(eigen)
pc2 <- eigen[2]/sum(eigen)

# Identify outliers - samples 4SDs from the mean
outliers <- as.character(plotdata$SampleID[c(which(plotdata$PC1 < mean(plotdata$PC1) - 4*sd(plotdata$PC1)),
                              which(plotdata$PC1 > mean(plotdata$PC1) + 4*sd(plotdata$PC1))), drop = T])

outliers <- c(outliers, as.character(plotdata$SampleID[c(which(plotdata$PC2 < mean(plotdata$PC2) - 4*sd(plotdata$PC2)),
                                           which(plotdata$PC2 > mean(plotdata$PC2) + 4*sd(plotdata$PC2))), drop = T] ))
  
plotdata <- left_join(plotdata, rownameToFirstColumn(covariates, "SampleID")) %>%
  dplyr::mutate(label = SampleID) %>% 
  dplyr::mutate(label = ifelse((label %in% outliers), label, NA))

# Bin Age of Death
plotdata <- plotdata %>%
  mutate(`AODbin` = ifelse(`Age of Death` < 60 & `Age of Death` >= 30, 60 , NA)) %>% 
  mutate(`AODbin` = ifelse(`Age of Death` < 30, 30, `AODbin`)) %>% 
  mutate(`AODbin` = ifelse(`Age of Death` < 90 & `Age of Death` >= 60, 89, `AODbin`)) %>%
  mutate(`AODbin` = ifelse(`Age of Death` == 90, 90, `AODbin`)) %>% 
  mutate(`AODbin` = as.character(`AODbin`)) %>% 
  mutate(`AODbin` = ifelse (`AODbin` == 30, "< 30", `AODbin`)) %>% 
  mutate(`AODbin` = ifelse (`AODbin` == 60, "< 60", `AODbin`)) %>% 
  mutate(`AODbin` = ifelse (`AODbin` == 89, "< 90", `AODbin`)) %>%
  mutate(`AODbin` = ifelse (`AODbin` == 90, "> 90", `AODbin`))

# Figure 2
p = list()
p[[1]] = ggplot(plotdata, aes(x = PC1, y = PC2, label = label))
p[[1]] = p[[1]] + geom_point(aes(color = Institution, shape = Dx, size = AODbin)) +
  ggtitle("PCA log2 CPM") + xlab("PC1 (89.9%)") + ylab("PC2 (4.6%)") + 
  theme_bw(12) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + 
  labs(shape = "Diagnosis", size = "Age of Death", tag = "a") +
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2), size = guide_legend(ncol = 2)) + 
  ggrepel::geom_text_repel(size = 3)

save_plot("./files/figure2a.pdf", p[[1]])

# Figure 2b
p[[2]] = ggplot(filt, aes (x= XIST, y = UTY)) 
p[[2]] = p[[2]] + geom_point(aes(color=`Sex`, shape = `Institution`)) + 
  ggtitle("Sex Check") + 
  theme_bw(12) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Sex", tag = "b")

save_plot("./files/figure2b.pdf", p[[2]], 
           base_aspect_ratio = 1)

