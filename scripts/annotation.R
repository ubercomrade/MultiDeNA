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
library(ReactomePA)
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
dirToWriteGO <- opt[["output_dir"]]

files <-  strsplit(files, ";")[[1]]
tools <-  strsplit(tools, ";")[[1]]
names(files) <- tools

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
  df <- df[df$annotation == "Promoter (1-2kb)" | df$annotation == "Promoter (1-2kb)" | df$annotation == "Promoter (2-3kb)",]
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
  for (i in tools) {
    df <- as.data.frame(peakAnnoList[i])
    names(df) <- colNames
    write.table(x = df,
                file = paste(writeDirectory, paste0(tolower(i),'_annotaion.tsv'), sep = '/'),
                sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
    
  }
}

dirToWriteReults <- '~/Downloads/'
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                      tssRegion=c(-3000, 3000), verbose=FALSE)
writeAnnotaion(peakAnnoList, dirToWriteReults)

genes <-  lapply(peakAnnoList, getGeneNamesByPromoters)
names(genes) <-  sub("_", "\n", names(genes))

pdf(paste(writeDirectory, 'tools_annotaion.pdf', sep = '/'))
plotAnnoBar(peakAnnoList)
dev.off()

enrich_go <- tryCatch(compareCluster(geneCluster  = genes,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           OrgDb = OrgDb),
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
if (is.S4(enrich_go)) {
  write.table(x = enrich_go@compareClusterResult,
              file = paste(writeDirectory, 'model_compare_GO.tsv', sep = '/'),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
  pdf(paste(writeDirectory, 'model_compare_GO.pdf', sep = '/'), width=10)
  dotplot(enrich_go, showCategory = 15, title = "Enrichment Analysis (GO)")
  dev.off()
}
#enrich_kegg <- tryCatch(compareCluster(geneCluster   = genes,
#                         fun           = "enrichKEGG",
#                         pvalueCutoff  = 0.05,
#                         pAdjustMethod = "BH"),
#                        error=function(cond) {
#                          message("There is no enrichment for current data")
#                          message("Here's the original error message:")
#                          message(cond)
#                          return(NA)
#                        },
#                        finally=function(cond) {
#                          return(NA)
#                        }
#                        )
#dotplot(enrich_kegg, showCategory = 15, title = "Enrichment Analysis (KEGG)")

enrich_pwm <- tryCatch(enrichGO(gene = genes$PWM,
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       ont = "BP",
                       OrgDb = OrgDb),
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

if (is.S4(enrich_pwm)) {
  write.table(x = enrich_pwm@result,
              file = paste(writeDirectory, 'pwm_GO.tsv', sep = '/'),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
  pdf(paste(writeDirectory, 'pwm_GO.pdf', sep = '/'), width=10)
  dotplot(enrich_pwm, showCategory = 20, title = "Enrichment Analysis (GO) for PWM sites")
  dev.off()
}

enrich_all <- tryCatch(enrichGO(gene = Reduce(c,genes),
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         ont = "BP",
                         OrgDb = OrgDb),
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

if (is.S4(enrich_all)) {
  write.table(x = enrich_pwm@result,
              file = paste(writeDirectory, 'all_GO.tsv', sep = '/'),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
  pdf(paste(writeDirectory, 'all_model_GO.pdf', sep = '/'), width=10)
  dotplot(enrich_all, showCategory = 20, title = "Enrichment Analysis (GO) for all models sites")
  dev.off()
}

