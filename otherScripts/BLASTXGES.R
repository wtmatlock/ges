library(readr)
library(dplyr)

df <- read_delim("blastx-output.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_names = FALSE, comment = "#", trim_ws = TRUE)

colnames(df) <- c('query.acc', 'subject.acc', 'perc.identity', 'alignment.length', 
                  'mismatches', 'gap.opens', 'q.start', 'q.end', 's.start', 
                  's.end', 'evalue', 'bit.score')

df <- distinct(df)
df <- df[df$`perc.identity`==100 & df$alignment.length==287,]
mult <- names(which(table(df$`query.acc`)>1))
df <- df[!(df$`query.acc` %in% mult),]

write.csv(paste0(df$query.acc, ".fasta"), '~/Desktop/contigs-GES-5.txt', col.names = FALSE, quote = FALSE, row.names = FALSE)