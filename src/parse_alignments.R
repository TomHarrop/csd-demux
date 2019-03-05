library(data.table)
library(ggplot2)

my_cn <- c("x", "subject_name", "subject_start", "subject_alnSize", "subject_strand", 
           "subject_seqSize", "query_name", "query_start", "query_alnSize", 
           "query_strand", "query_seqSize", "alignment", "score1", "score2")
           

aln <- fread("grep -v '^#' output/020_mapped/aln.tab", col.names = my_cn)

aln[, length(unique(query_name)), by = subject_name]
aln[, length(unique(query_name))]

pd1 <- unique(aln, by = c("query_name", "subject_name"))
ggplot(pd1, aes(x = query_seqSize)) +
    xlim(c(100, 1500)) +
    facet_wrap(~ subject_name) +
    geom_histogram(binwidth = 10)
