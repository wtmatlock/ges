# mash sketch -s 1000000 ./contigs-blastx/*.fasta -o all-blastx
# echo contigs-blastx/*.fasta | xargs -n 1 mash dist all-blastx.msh > mash_edgelist.tsv

library(readr)

# n x n distance matrix as edgelist
elist <- read_delim("~/Desktop/GES/data/mash_edgelist.tsv", 
                                     delim = "\t", escape_double = FALSE, 
                                     col_names = FALSE, trim_ws = TRUE)
colnames(elist) <- c("contig1", "contig2", "mdist", "pval", "jaccard")

# format contig names
elist$contig1 <- gsub("./contigs-GES-5/", "", elist$contig1)
elist$contig2 <- gsub("contigs-GES-5/", "", elist$contig2)

elist2 <- elist
original <- unique(elist$contig1) # n=431

# filter edgelist for duplicates
elist <- elist[elist$mdist<0.0001,]
elist<- elist[elist$contig1!=elist$contig2,]
elist <- elist[1:3]

# bioprojects

NCBI_GES <- read_delim("~/Desktop/GES/data/NCBI-GES.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

write.csv(NCBI_GES$Contig, '~/Desktop/NCBI-GES-contig.txt', row.names = FALSE, quote=FALSE)

# select representative from connected components 

library(igraph)
library(dplyr)

g <- graph.data.frame(d=elist, directed=FALSE)
plot(g, vertex.label = V(g)$name)
g_cc <- as.data.frame(components(g)[1])

duplicates <- rownames(g_cc) # n=258
g_cc <- g_cc %>%
  add_rownames("contig") %>%
  group_by(membership) %>%
  mutate(rep = first(contig))

keep <- unique(g_cc$rep) # n=56
remove <- setdiff(duplicates, keep) # n=202
new <- setdiff(original, remove)

# check 

elist2 <- elist2[elist2$contig1 %in% new,]
elist2 <- elist2[elist2$contig2 %in% new,]
elist2 <- elist2[elist2$mdist<0.0001,]
elist2<- elist2[elist2$contig1!=elist2$contig2,] # should be empty!

write.csv(new, '~/Desktop/contigs-GES-5-dedup.txt', col.names = FALSE, quote = FALSE, row.names = FALSE)