# Viz HiChip loops

# Code
library(tidyverse)
library(plotgardenerData)

# Paths
path_main <- '/Users/IEO5505/Desktop/fragile_enhancer_clinical'
path_data <- paste0(path_main, '/data/Hi_Chip/loops/')
path_results <- paste0(path_main, '/results/Hi_Chip/')

# Read data
SCR_loops <- read.table(paste0(path_data, '/SCR_loops.bed'), header=TRUE)
SCR_loops$loop <- paste(SCR_loops$chr1,SCR_loops$s1, SCR_loops$e1, 'chr', SCR_loops$chr2, SCR_loops$s2, SCR_loops$e2, sep='_')
SCR_loops$loop <- paste('chr', SCR_loops$loop, sep='')
enh_loops <- read.csv(paste0(path_results, '/SCR_CtIP_enh_loops.csv'), row.names = 1)

# Viz some loops
enh_loops %>% arrange(Q.Value_Bias)

# Prova
chrom = "chr20"
chromstart = 348000
chromend = 395000

# Plot loop
plotBedpe(SCR_loops, chrom, chromstart, chromend, heights = SCR_loops$Dist, plottype="loops", flip=TRUE)
labelgenome(chrom, chromstart,chromend,side=3, n=3,scale="Mb")
axis(side=2,las=2,tcl=.2)
mtext("contact freq",side=2,line=1.75,cex=.75,font=2)

