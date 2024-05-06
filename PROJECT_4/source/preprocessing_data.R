pacman::p_load(dplyr,DESeq2)
path <- "/data/gpfs-1/users/phgi10_c/work/Uni/AML/PROJECT4/data/GSE68086_TEP_data_matrix.txt"
raw_data <- read.csv(path,sep="\t",header=TRUE) %>% as_tibble()

filter_data <- raw_data %>% 
  filter(across(everything(), ~ . > 5))

#Create count matrix
count_matrix_aml <- filter_data[,1:ncol(filter_data)]
count_mt <- count_matrix_aml[,2:ncol(count_matrix_aml)]
#rownames(count_matrix_aml) <- filter_data$X
#Create metadata
pattern <- c("HD", "control")
samples <- colnames(filter_data)
samples <- samples[2:length(samples)]
dex <- ifelse(grepl("HD|Control", samples), "control", "treatment")
metadata <- data.frame(id = samples,dex = dex)

#DESeq2
dds <- DESeqDataSetFromMatrix(countData=count_mt,
                              colData=metadata,
                              design=~dex)

#####################
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
#normalizedCounts <- t(t(counts(dds)) / dds$sizeFactors)
#keep <- filterByExpr(dds)
#dds <- dds[keep,]\
res <- results(dds)
#res <- nbinomWaldTest(dds, "control", "treatment")
res["genes"] <- filter_data$X
sig_indices <-which(res$padj < 0.05 & (res$log2FoldChange < -1 | res$log2FoldChange > 1))

normalizedCounts <- counts(dds, normalized = TRUE)
df <- data.frame(genes=filter_data$X,normalizedCounts,padj=res$padj,log2fc=res$log2FoldChange)
significantGenes <- df[sig_indices,]
significantGenes$genes
dim(significantGenes)
significantGenes <- significantGenes[, -which(names(significantGenes) == "padj" | names(significantGenes) == "log2fc")]
write.csv(significantGenes, "/data/gpfs-1/users/phgi10_c/work/Uni/AML/PROJECT4/preprocessed_data.csv", row.names = FALSE)
