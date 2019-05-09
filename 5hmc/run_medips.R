library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)

my_args <- commandArgs(T)
bam_path <- my_args[1]
wig_path <- my_args[2]
sample <- my_args[3]

MSet <- MEDIPS.createSet(file=bam_path, extend=0, shift=0, window_size=1000, BSgenome="BSgenome.Hsapiens.UCSC.hg19", uniq=1e-3)
MEDIPS.exportWIG(Set=MSet, file=wig_path, format="count", descr=sample)

