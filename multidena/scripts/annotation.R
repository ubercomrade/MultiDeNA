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
#library(ChIPseeker)
library(clusterProfiler)
library(optparse)

option_list = list(
  make_option(c("--input_annotations"), dest="input_annotations", action="store",
              help="separated list of files (sep = ;)"),
  make_option(c("--models_names"), dest="models_names", action="store",
              help="separated list of names (sep = ;)"),
  make_option(c("--genome"), dest="genome", action="store",
              help="genome version from list [mm10, hg38, tair10]"),
  make_option(c("--output_dir"), dest="output_dir", action="store",
              help="directory name to write results"))

opt = parse_args(OptionParser(option_list=option_list))

files <-  opt[["input_annotations"]]
tools <-  opt[["models_names"]]
genome <-  opt[["genome"]]
writeDirectory <- opt[["output_dir"]]

#files <- "/home/anton/Documents/PhD/gtrd-mm10/results/PEAKS059242_MYOD1_P10085_MACS2/annotation_vg/test/pwm_test_1.00e-04.ann_genes.txt;/home/anton/Documents/PhD/gtrd-mm10/results/PEAKS059242_MYOD1_P10085_MACS2/annotation_vg/test/bamm_test_1.00e-04.ann_genes.txt;/home/anton/Documents/PhD/gtrd-mm10/results/PEAKS059242_MYOD1_P10085_MACS2/annotation_vg/test/sitega_test_1.00e-04.ann_genes.txt"
#tools <- "PWM;BaMM;SiteGA"
#genome <- "mm10"
#writeDirectory <- '/home/anton/Downloads/test/'

files <-  strsplit(files, ";")[[1]]
tools <-  strsplit(tools, ";")[[1]]
names(files) <- tools

if (genome == "tair10"){
  keyType = "TAIR"
} else {
  keyType = "ENSEMBL"
}


if (genome == 'hg38') {
  OrgDb <- 'org.Hs.eg.db'
}

if (genome == 'mm10') {
  OrgDb <- 'org.Mm.eg.db'
}

if (genome == 'tair10') {
  OrgDb <- 'org.At.tair.db'
}

readGenes <- function(path){
  return(scan(path, what="", sep="\n"))
}
genes <-  lapply(files, readGenes)

print(genes)
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