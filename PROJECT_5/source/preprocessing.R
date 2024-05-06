library(readr)
library(GEOquery)
library(dplyr)
require(hta20transcriptcluster.db)
#read the top significant genes found by using GEO2R tool:
GSE100179_top_table <- read_delim("GSE100179.top.table.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
#filter with the help of padj. values
resFilt_GSE100179 <- GSE100179_top_table[which(GSE100179_top_table$adj.P.Val < 0.05), ]


# load series and platform data from GEO
gset1 <- getGEO("GSE100179", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
# choose the index
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]
#load the expression data.
t0<- gset1@assayData[["exprs"]]

# extract the information of cancer/normal/adenoma from the metadata.
title0<- gset1@phenoData@data[["characteristics_ch1"]]
#extract the sample names from the expression data.
clname<- colnames(t0)
# make dataframe of the sample names and their corresponding information about cancer/non cancer/adema
d_col<- as.data.frame(title0,clname)

gene_data_frame = fData(gset1)

#map the IDs of expression data and convert them to gene symbols:
mapping <- select(
  hta20transcriptcluster.db,
  keys = rownames(gset1),
  column = c('SYMBOL', 'ENTREZID', 'ENSEMBL'),
  keytype = 'PROBEID')

mapped<- (mapping[!is.na(mapping$SYMBOL),])
#merge the mapped genes symbols with expression data
d<-merge(t0,mapped, by.x=0, by.y= "PROBEID")
d_f<- merge(d,resFilt_GSE100179, by.x="Row.names", by.y= "ID")
# make gene symbols unique:
d_f1<- d_f %>% distinct(SYMBOL, .keep_all = T)

#################do the same for another data##############################################
GSE117606_top_table <- read_delim("GSE117606.top.table.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
resFilt_GSE117606 <- GSE117606_top_table[which(GSE117606_top_table$adj.P.Val < 0.05), ]

gset1 <- getGEO("GSE117606", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]
t0<- gset1@assayData[["exprs"]]

gene_data_frame = fData(gset1)


title0<- gset1@phenoData@data[["title"]]
clname<- colnames(t0)

a_col<- as.data.frame(title0,clname)

library(data.table)
names<- fread("names.txt")
a<-merge(t0,names, by.x=0, by.y= "ID")

a_f<- merge(a,resFilt_GSE117606, by.x="Row.names", by.y= "ID")
library(dplyr)
a_f1<- a_f %>% distinct(GeneSymbol, .keep_all = T)

###################################################################
GSE4183_top_table <- read_delim("GSE4183.top.table.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

resFilt_GSE4183 <- GSE4183_top_table[which(GSE4183_top_table$adj.P.Val < 0.05), ]#


gset1 <- getGEO("GSE4183", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]
t0<- gset1@assayData[["exprs"]]

gene_data_frame = fData(gset1)


title0<- gset1@phenoData@data[["title"]]
clname<- colnames(t0)

b_col<- as.data.frame(title0,clname)

b<-merge(t0,gene_data_frame, by.x=0, by.y= "ID")

b_f<- merge(b,resFilt_GSE4183, by.x="Row.names", by.y= "ID")
b_f1<- b_f %>% distinct(`Gene symbol`, .keep_all = T)

#############################################################
GSE71187_top_table <- read_delim("GSE71187.top.table.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
resFilt_GSE71187 <- GSE71187_top_table[which(GSE71187_top_table$adj.P.Val < 0.05), ]#

gset1 <- getGEO("GSE71187", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]
t0<- gset1@assayData[["exprs"]]


title0<- gset1@phenoData@data[["title"]]
clname<- colnames(t0)

cc_col<- as.data.frame(title0,clname)


gene_data_frame = fData(gset1)

cc<-merge(t0,gene_data_frame, by.x=0, by.y= "ID")

cc_f<- merge(cc,resFilt_GSE71187, by.x="Row.names", by.y= "ID")
cc_f1<- cc_f %>% distinct(Gene.symbol, .keep_all = T)
###############################################################################
#merge all the data with the help of genesymbols.
ab<- merge(a_f1,b_f1,by.x="GeneSymbol",by.y="Gene symbol")
abc<- merge(ab,cc_f1,by.x = "GeneSymbol", by.y= "Gene symbol")
abcd <- merge(abc,d_f1,by.x="GeneSymbol",by.y= "SYMBOL")

############################################################
#order the genesymbol with help of a column:
occurrence <- abcd[order(abcd$GeneSymbol, abcd$GSM3304749, decreasing=TRUE),]
#clean any duplicate gene symbols:
occurrenceClean <- occurrence[!duplicated(occurrence$GeneSymbol),]
cleany<- occurrenceClean[ , colSums(is.na(occurrenceClean))==0]# remove columns with NA
# remove last 5 unnecessary columns
cleany<- cleany[1:(length(cleany)-5)]
#again remove more columns
last_cleany<- cleany[,-2]
last_cleany<- na.omit(last_cleany)# remove any rows with NA.
#bind all metadata cols
binded_col<- rbind(a_col,b_col,cc_col,d_col)
#write it and save:
fwrite(binded_col, file = "metadata_col_info.csv", sep = ",",row.names = TRUE)
fwrite(last_cleany, file = "cleanedpreproc_data.csv",sep= ",")

