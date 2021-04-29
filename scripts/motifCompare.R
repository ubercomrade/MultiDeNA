#!/usr/bin/env Rscript
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("motifStack")
library(motifStack)

option_list = list(
  make_option(c("--dir_with_motifs"), dest="dir_with_motifs", action="store",
              help="dir with meme motifs"),
  make_option(c("--dir_to_write"), dest="dir_to_write", action="store",
            help="dir to write results"))
opt = parse_args(OptionParser(option_list=option_list))

directory <-  opt[["dir_with_motifs"]]
writeDirectory <- opt[["dir_to_write"]]

files <- dir(file.path(directory),"meme$", full.names = TRUE)
files <- files[!grepl('sitega', files)]
motifs <- importMatrix(files, format='meme')

tools <- c('PWM', 'diPWM', 'BaMM', 'InMoDe', 'StruM')
names(tools) <- list('pwm', 'dipwm', 'bamm', 'inmode', 'strum')

for (i in names(motifs)) {
  motifs[[i]]@name = tools[[i]]
}

pdf(paste(writeDirectory, 'compare_motif_logo.pdf', sep = '/'), width=10)
motifStack(motifs, layout="tree")
dev.off()
