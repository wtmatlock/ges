require(ggplot2)
require(cowplot)
require(gggenes)
require(ggdendro)
require(reshape2)
require(vegan) # for Jaccard distance

arg1 <- '/Users/willmatlock/Desktop/GES/data/pangraph-output-5kbp/pangraph-export.gfa.blocks.csv'
arg2 <- 'NLGAUUIIFT'

# Genome blocks approach
genome.blocks <- read.csv(arg1, header=T, stringsAsFactors = F)
genome.blocks <- genome.blocks[!(genome.blocks$genome %in% c('NZ_VAAM01000031.1', 'NZ_CP073313.1')),] # hacky bit
#colnames(genome.blocks) <- c("genome", "block", "strand", "start", "end", "colour")
genome.blocks$forward <- ifelse(genome.blocks$strand=="+", TRUE, FALSE)
block.counts <- table(genome.blocks$block)
blocks.which.need.colours <- names(block.counts)[which(block.counts>1)]

genome.blocks$block.coloured <- sapply(genome.blocks$block,
                                       function(x) ifelse(x %in% blocks.which.need.colours,
                                                          x,
                                                          "_other"))
block.colours <- unique(genome.blocks$colour)
names(block.colours) <- unique(genome.blocks$block.coloured)

GES.block.locations <- genome.blocks[which(genome.blocks$block==arg2),c("start", "end")]
rownames(GES.block.locations) <- genome.blocks[which(genome.blocks$block==arg2),"genome"]
transformAnnotationsBlocks <- function(genome, default_offset=10000){
  offset <- default_offset-GES.block.locations[genome, c("start")]
  return(offset)
}

genome.blocks$new.start <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["start"])+transformAnnotationsBlocks(x["genome"]))-10001
genome.blocks$new.end <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["end"])+transformAnnotationsBlocks(x["genome"]))-10001

# Get all the genome paths (in terms of blocks)
genome.paths <- sapply(unique(genome.blocks$genome), function(x) paste(genome.blocks[which(genome.blocks$genome==x), "block"], collapse=","))
genome.paths <- sort(genome.paths, decreasing = TRUE)

# genomes with just gene block
genome.paths[genome.paths =="NLGAUUIIFT"]

# Subset to unique paths
# Keep one representative genome for each
genome.path.reps <- names(genome.paths)[!duplicated(genome.paths)]
genome.blocks.unique <- genome.blocks[which(genome.blocks$genome %in% genome.path.reps),]
# And store the number of examples of each
genome.blocks.unique$genome.path <- genome.paths[genome.blocks.unique$genome]
genome.blocks.unique$n.reps <- sapply(genome.blocks.unique$genome.path,
                                    function(x) table(genome.paths)[x])
# All genomes for each rep
genome.blocks.unique$all.reps <- sapply(genome.blocks.unique$genome.path,
                                      function(x) names(which(genome.paths==x)))
genome.blocks.unique$genome.path.name <- paste0("Type", as.numeric(as.factor(genome.blocks.unique$genome.path )))
genome.blocks.unique$genome.n <- sapply(genome.blocks.unique$n.reps, function(x) ifelse(x==1,
                                                                                        "", paste0("n=", x)))

# Now use a dendrogram to order the genomes
m <- acast(genome ~ block, data=genome.blocks.unique, fill=0, fun.aggregate=length)
m.dist <- vegdist(m, method="jaccard") # jaccard distances based on block presence/absence
dendro <- as.dendrogram(hclust(m.dist))
dendro_order <- order.dendrogram(dendro)
genome_labels <- dendro_data(dendro)$labels$label
genome.blocks.unique$genome.ordered <- ordered(genome.blocks.unique$genome,
                                               levels=genome_labels)
# Position from bottom to top
position <- cbind(genome_labels, length(genome_labels):1)
plot.placement <- merge(genome.blocks.unique, position, by.x='genome.ordered', by.y='genome_labels', all.x=TRUE)
plot.placement <- plot.placement[,c("genome.ordered", "all.reps", "V2")]
colnames(plot.placement) <- c("genome.ordered", "all.reps", "pos")
library(tidyr)
plot.placement <- unnest(plot.placement, cols=all.reps)
plot.placement <- unique(plot.placement)

# Make the linear block plot
p.blocks <- ggplot(genome.blocks.unique, aes(xmin = new.start, xmax = new.end, forward = forward, y = genome.ordered, fill = block.coloured)) +
  geom_gene_arrow(arrow_body_height = unit(1, "mm"),
                  arrowhead_height = unit(1, "mm"),
                  arrowhead_width = unit(1, "mm")) +
  theme_genes()+
  scale_fill_manual(values=block.colours)+
  ylab("")+
  theme(legend.position = "none")+
  scale_y_discrete(breaks=genome.blocks.unique$genome.ordered, labels=genome.blocks.unique$genome.n)+
  ggtitle("Linearised blocks")+
  theme(plot.title=element_text(hjust=0.5))+
  xlab("Position (bp)")

p.blocks

# import duplication data

# wrangling
library(readr)
dup <- read_csv("Desktop/GES/data/rmdup.txt", col_names = FALSE)
library(stringr)
dup$original.size <- word(dup$X1, 1, sep="\t")
dup$X1 <- word(dup$X1, 2, sep="\t")
dup$group <- 1:nrow(dup)
dup <- dup[,c(1:27, 36:37)]
dup <- melt(dup, id=c("original.size","group"), value.name='original.genome')
dup <- drop_na(dup)
dup <- dup[c("group", "original.genome", "original.size")]

library(dplyr)
library(stringr)

genome.blocks.unique.melt <- unnest(genome.blocks.unique, all.reps)
genome.blocks.unique.melt <- merge(genome.blocks.unique.melt, dup, by.x='all.reps', by.y='original.genome', all.x=TRUE)
genome.blocks.unique.melt <- genome.blocks.unique.melt[c('genome.ordered', 'genome.path.name', 'n.reps', 'original.size')]
genome.blocks.unique.melt <- unique(genome.blocks.unique.melt)
genome.blocks.unique.melt$original.size <- as.numeric(genome.blocks.unique.melt$original.size)
genome.blocks.unique.melt$genome.path.name <- as.factor(genome.blocks.unique.melt$genome.path.name)

genome.blocks.unique.melt[genome.blocks.unique.melt$original.size>0&!(is.na(genome.blocks.unique.melt$original.size)),]$original.size <- 
  genome.blocks.unique.melt[genome.blocks.unique.melt$original.size>0&!(is.na(genome.blocks.unique.melt$original.size)),]$original.size -1

genome.blocks.unique.melt <- genome.blocks.unique.melt %>%
  group_by(genome.path.name) %>%
  mutate(original.total = sum(original.size, na.rm = TRUE))

removed <- genome.blocks.unique.melt[c("genome.ordered", "n.reps", "original.total")]
colnames(removed) <- c("genome.ordered", "kept", 'removed')
removed <- unique(removed)

removed <- merge(removed, plot.placement[c("genome.ordered", 'pos')], all.x=TRUE)
removed <- unique(removed)
removed$pos <- as.numeric(removed$pos)
removed$plot.pos <- abs(removed$pos - max(removed$pos)) +1

sum(removed$kept)
431-5-sum(removed$removed)

removed <- melt(removed, id.vars = c('genome.ordered', 'pos', 'plot.pos'))
removed$variable <- factor(removed$variable, levels = c("removed", "kept"))

p.rm <- ggplot(removed, aes(x=plot.pos, y=value)) + 
  geom_col(aes(y=value, fill=variable),size=0.25, color='black', width = 1) + 
  coord_flip() +
  scale_x_discrete(labels = NULL, breaks = removed$plot.pos) +
  labs(x = "", y='Count') +
  scale_fill_manual("", labels=c("Removed", "Kept"), values = c("#e41a1c","#377eb8")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Deduplication") #+ theme(legend.position = "none")



library(cowplot)
grid::current.viewport()
cowplot::plot_grid(p.blocks, p.rm,
                   nrow=1, align='h', rel_widths = c(1, 1),
                   labels=c("a", "b"))



