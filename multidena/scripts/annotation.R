#!/usr/bin/env Rscript
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ChIPseeker")
#BiocManager::install("clusterProfiler")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("org.At.tair.db")
#BiocManager::install("ReactomePA")

#remotes::install_version("RSQLite", version = "2.2.5")
library(ChIPseeker)
library(clusterProfiler)
library(optparse)

option_list = list(
  make_option(c("--input_scans"), dest="input_scans", action="store",
              help="separated list of files (sep = ;)"),
  make_option(c("--models_names"), dest="models_names", action="store",
              help="separated list of names (sep = ;)"),
  make_option(c("--genome"), dest="genome", action="store",
              help="genome version from list [mm10, hg38, tair10]"),
  make_option(c("--output_dir"), dest="output_dir", action="store",
              help="directory name to write results"))

opt = parse_args(OptionParser(option_list=option_list))

files <-  opt[["input_scans"]]
tools <-  opt[["models_names"]]
genome <-  opt[["genome"]]
writeDirectory <- opt[["output_dir"]]

#files <- "/Users/anton/Documents/PhD/gtrd-tair10-choosen/results_sitega//PEAKS042881_CCA1_P92973_MACS2_1344/scan//pwm_test_1.90e-04.bed;/Users/anton/Documents/PhD/gtrd-tair10-choosen/results_sitega//PEAKS042881_CCA1_P92973_MACS2_1344/scan//bamm_test_1.90e-04.bed;/Users/anton/Documents/PhD/gtrd-tair10-choosen/results_sitega//PEAKS042881_CCA1_P92973_MACS2_1344/scan//inmode_test_1.90e-04.bed;/Users/anton/Documents/PhD/gtrd-tair10-choosen/results_sitega//PEAKS042881_CCA1_P92973_MACS2_1344/scan//sitega_test_1.90e-04.bed"
#tools <- "PWM;BaMM;InMoDe;SiteGA"
#genome <- "tair10"

files <-  strsplit(files, ";")[[1]]
tools <-  strsplit(tools, ";")[[1]]
names(files) <- tools

if (genome == "tair10"){
  keyType = "TAIR"
} else {
  keyType = "ENTREZID"
}


if (genome == 'hg38') {
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  OrgDb <- 'org.Hs.eg.db'
}

if (genome == 'mm10') {
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  OrgDb <- 'org.Mm.eg.db'
}

if (genome == 'tair10') {
  library(TxDb.Athaliana.BioMart.plantsmart28)
  library(org.At.tair.db)
  txdb <- TxDb.Athaliana.BioMart.plantsmart28
  OrgDb <- 'org.At.tair.db'
}

getGeneNamesByPromoters <- function(i) {
  df <- as.data.frame(i)
  df <- df[df$annotation == "Promoter" | df$annotation == "Promoter (1-2kb)" | df$annotation == "Promoter (2-3kb)",]
  genes <- unique(df$geneId)
  return(genes)
}

writeAnnotaion <- function(peakAnnoList, writeDirectory) {
  tools <-  names(peakAnnoList)
  colNames <- c("chromosome", "start", "end", "width",
                "str", "name", "score", "strand",
                "site","annotation", "geneChr",
                "geneStart", "geneEnd", "geneLength",
                "geneStrand", "geneId", "transcriptId",
                "distanceToTSS")
  colNames <- c("chromosome", "start",  "end", "width",
                "strand", "name", "score", "site", "annotation",
                "geneChr", "geneStart",  "geneEnd",
                "geneLength", "geneStrand", "geneId",
                "transcriptId", "distanceToTSS")
  for (i in tools) {
    df <- as.data.frame(peakAnnoList[i])
    names(df) <- colNames
    write.table(x = df,
                file = paste(writeDirectory, paste0(tolower(i),'_annotaion.tsv'), sep = '/'),
                sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )

  }
}

CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

readScanTablesTair <- function(path){
  df <- read.table(path)
  colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "site")
  #df$chr <- sapply(df$chr, CapStr)
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  seqlevelsStyle(gr) <- "Ensembl"
  return(gr)
}

readScanTablesMammals <- function(path){
  df <- read.table(path)
  colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "site")
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  seqlevelsStyle(gr) <- "UCSC"
  return(gr)
}

if (genome == "tair10"){
  grs <- lapply(files, readScanTablesTair)
  peakAnnoList <- lapply(grs, annotatePeak, TxDb=txdb,
                         tssRegion=c(-1000, 1000), verbose=FALSE)

} else {
  grs <- lapply(files, readScanTablesMammals)
  peakAnnoList <- lapply(grs, annotatePeak, TxDb=txdb,
                         tssRegion=c(-3000, 3000), verbose=FALSE)
}

writeAnnotaion(peakAnnoList, writeDirectory)
genes <-  lapply(peakAnnoList, getGeneNamesByPromoters)
#genes <-  lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
#genes <-  lapply(grs, function(i) seq2gene(i, tssRegion=c(-1000, 1000), flankDistance = 3000, TxDb=txdb))
names(genes) <-  sub("_", "\n", names(genes))

pdf(paste(writeDirectory, 'tools_annotaion.pdf', sep = '/'))
plotAnnoBar(peakAnnoList)
dev.off()


enrich_go <- tryCatch(compareCluster(geneCluster  = genes,
                                     fun           = "enrichGO",
                                     pvalueCutoff  = 0.05,
                                     qvalueCutoff = 0.05,
                                     ont = "BP",
                                     pAdjustMethod = "BH",
                                     OrgDb = OrgDb,
                                     keyType = keyType),
                      error=function(cond) {
                        message("There is no enrichment for current data")
                        message("Here's the original error message:")
                        message(cond)
                        return(NA)
                      },
                      finally=function(cond) {
                        return(NA)
                      }
)

#l <- length(unique(as.data.frame(enrich_go)$ID))
#dotplot(enrich_go, showCategory = 42, title = "Enrichment Analysis (GO)", font.size = 10, label_format = 70)

if (isS4(enrich_go)) {
  df <- enrich_go@compareClusterResult
  df <- df[df$p.adjust <= 0.05,]
  write.table(x = df,
              file = paste(writeDirectory, 'compare_models_GO.tsv', sep = '/'),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
  l <- length(unique(as.data.frame(enrich_go)$ID))
  p <- dotplot(enrich_go, showCategory = l, title = "Enrichment Analysis (GO)", font.size = 10, label_format = 50)
  pdf(paste(writeDirectory, 'compare_models_GO.pdf', sep = '/'), width=10)
  print(p)
  dev.off()
} else {
  colNames <- c("Cluster", "ID", "Description", "GeneRatio",
                "BgRatio", "pvalue", "p.adjust", "qvalue",
                "geneID", "Count")
  df <- data.frame(matrix(ncol=10,nrow=0,
                          dimnames=list(NULL, colNames)))
  write.table(x = df,
              file = paste(writeDirectory, 'compare_models_GO.tsv', sep = '/'),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
}

enrich_pwm <- tryCatch(enrichGO(gene = genes$PWM,
                                pvalueCutoff  = 0.05,
                                qvalueCutoff = 0.05,
                                ont = "BP",
                                pAdjustMethod = "BH",
                                OrgDb = OrgDb,
                                keyType = keyType),
                       error=function(cond) {
                         message("There is no enrichment for current data")
                         message("Here's the original error message:")
                         message(cond)
                         return(NA)
                       },
                       finally=function(cond) {
                         return(NA)
                       }
)

#dotplot(enrich_pwm, showCategory = 20, title = "Enrichment Analysis (GO) for PWM sites")

if (isS4(enrich_pwm)) {
  df <- enrich_pwm@result
  df <- df[df$p.adjust <= 0.05,]
  l <- length(df$ID)
  write.table(x = df,
              file = paste(writeDirectory, 'pwm_model_GO.tsv', sep = '/'),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
  p <- dotplot(enrich_pwm, showCategory = l, title = "Enrichment Analysis (GO) for PWM sites")
  pdf(paste(writeDirectory, 'pwm_model_GO.pdf', sep = '/'), width=10)
  print(p)
  dev.off()
}

enrich_all <- tryCatch(enrichGO(gene = Reduce(c,genes),
                                pvalueCutoff  = 0.05,
                                qvalueCutoff = 0.05,
                                ont = "BP",
                                pAdjustMethod = "BH",
                                OrgDb = OrgDb,
                                keyType = keyType),
                       error=function(cond) {
                         message("There is no enrichment for current data")
                         message("Here's the original error message:")
                         message(cond)
                         return(NA)
                       },
                       finally=function(cond) {
                         return(NA)
                       }
)

#dotplot(enrich_all, showCategory = 20, title = "Enrichment Analysis (GO) for all models sites")

if (isS4(enrich_all)) {
  df <- enrich_all@result
  df <- df[df$p.adjust <= 0.05,]
  l <- length(df$ID)
  write.table(x = df,
              file = paste(writeDirectory, 'all_models_GO.tsv', sep = '/'),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
  p <- dotplot(enrich_all, showCategory = l, title = "Enrichment Analysis (GO) for all models sites")
  pdf(paste(writeDirectory, 'all_models_GO.pdf', sep = '/'), width=10)
  print(p)
  dev.off()
}
