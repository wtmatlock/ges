library(readr)
library(dplyr)

lengths <- read_delim("contigs-GES-5-lengths.tsv", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
colnames(lengths) <- c("contig", "length")

ges <- read_delim("blastx-output.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_names = FALSE, comment = "#", trim_ws = TRUE)
colnames(ges) <- c('contig', 'subject.acc', 'perc.identity', 'alignment.length', 
                  'mismatches', 'gap.opens', 'q.start', 'q.end', 's.start', 
                  's.end', 'evalue', 'bit.score')

ges <- distinct(ges)
ges <- ges[ges$'perc.identity'==100 & ges$alignment.length==287,]
mult <- names(which(table(ges$contig)>1))
ges <- ges[!(ges$contig %in% mult),]

ges <- ges[c("contig", "q.start", "q.end")]

df <- merge(ges, lengths, by.x='contig', by.y='contig', all.x=TRUE)

df$strand <- ifelse(df$q.start>df$q.end, "+", "-")

df <- df[!(df$query.acc %in% names(which(table(df$query.acc)>1))),]

write.csv(paste0(df$query.acc, ".fasta"), '~/Desktop/contigs-GES-5.txt', col.names = FALSE, quote = FALSE, row.names = FALSE)