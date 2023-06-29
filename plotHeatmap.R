library(readr)

snps <- read_delim("/Users/willmatlock/Desktop/integrase-snps/snp-matrix.tab", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
names <- gsub("-integrase", "", snps$`snp-dists 0.8.2`)
snps <- snps[,-1]
colnames(snps) <- names
rownames(snps) <- names
snps <- as.matrix(snps)


library("pheatmap")
plot <- pheatmap(snps,display_numbers = snps, cutree_rows = 2, cutree_cols = 2,
         clustering_distance_cols="euclidean", clustering_method="complete")

dendro <- sort(cutree(plot$tree_row, k=2))

write.csv(as.data.frame(dendro), '~/Desktop/integrase-clusters.csv')

c1 <- names(dendro[dendro==1])
length(c1)
c2 <- names(dendro[dendro==2])
length(c2)

snps.c1 <- snps[c1, c1]
snps.c1 <- snps.c1[upper.tri(snps.c1, diag=FALSE)]
summary(snps.c1)

snps.c2 <- snps[c2, c2]
snps.c2 <- snps.c2[upper.tri(snps.c2, diag=FALSE)]
summary(snps.c2)
