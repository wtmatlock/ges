require(readr)
require(ggplot2)
require(cowplot)
require(gggenes)
require(ggdendro)
require(reshape2)
require(vegan) # for Jaccard distance

arg1 <- '/Users/willmatlock/Desktop/pangraph-output-10000/pangraph-export.gfa.blocks.csv'
arg2 <- 'GWKFGRWJMP'

# Genome blocks approach
genome.blocks <- read.csv(arg1, header=T, stringsAsFactors = F)
genome.blocks <- genome.blocks[!(genome.blocks$genome %in% c('NZ_CP073313.1', 'NZ_VAAM01000031.1', 'NZ_UARQ01000013.1')),] # hacky bit
contigs_blastx_dedup <- read_csv("Desktop/contigs-blastx-dedup.txt", 
                                 col_names = FALSE)
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
genome.paths[genome.paths =="HDYRKRALIV"]

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
  geom_gene_arrow(arrow_body_height = unit(2, "mm"),
                  arrowhead_height = unit(2, "mm"),
                  arrowhead_width = unit(1, "mm")) +
  theme_genes()+
  scale_fill_manual(values=block.colours)+
  ylab("")+
  theme(legend.position = "none")+
  scale_y_discrete(breaks=genome.blocks.unique$genome.ordered, labels=genome.blocks.unique$genome.n)+
  ggtitle("Linearised blocks")+
  theme(plot.title=element_text(hjust=0.5))+
  xlab("Position (bp)")

#p.blocks

# plot dendrogram

p.dendro <- ggdendrogram(dendro) + 
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  theme_dendro()

#p.dendro

library(ggdendro)

# import metadata

library(dplyr)
library(stringr)

metadata <- read_delim("~/Desktop/GES/data/NCBI-GES.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

df <- merge(genome.blocks, metadata, by.x='genome', by.y='Contig', all.x=TRUE)
df <- df[,c('genome', '#Scientific name', 'Location', 'Collection date')]
df <- unique(df)
df <- merge(df, plot.placement, by.x='genome', by.y='all.reps', all.x=TRUE)
colnames(df) <- c("genome", "genus", "country", "sampling.year", "genome.ordered", "pos")
df$pos <- as.numeric(df$pos)
df$plot.pos <- abs(df$pos - max(df$pos)) +1


df$genus <- word(df$genus, 1)
df$country <- word(df$country, 1, sep=":")
df$sampling.year <- substr(df$sampling.year, 1, 4)

who_region <- read_csv("~/Desktop/GES/data/who-region.csv")
who_region <- who_region[,c("name", 'region')]
df <- merge(df, who_region, by.x='country', by.y='name', all.x=TRUE)
df$region <- ifelse(df$country=='Czech Republic', "Europe",
                    ifelse(df$country=="Russia", "Europe",
                           ifelse(df$country=="Taiwan", "Asia",
                                  ifelse(df$country=="United Kingdom", 'Europe',
                                         ifelse(df$country=="USA", "Americas", df$region))))) 


p.genus <- ggplot(df, aes(x=plot.pos)) + 
  geom_bar(aes(y=..count.., fill=genus), position='fill',size=0.25, color='black', width = 1) + 
  coord_flip() +
  scale_x_discrete(labels = NULL, breaks = df$plot.pos) +
  labs(x = "", y='Proportion') +
  scale_fill_brewer("", palette = "Set3") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Genus") + theme(legend.position = "none")

#p.genus

p.loc <- ggplot(df, aes(x=plot.pos)) + 
  geom_bar(aes(y=..count.., fill=region), position='fill',size=0.25, color='black', width = 1) + 
  coord_flip() +
  scale_x_discrete(labels = NULL, breaks = df$plot.pos) +
  labs(x = "", y='Proportion') +
  scale_fill_brewer("", palette="Accent") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("WHO region") + theme(legend.position = "none")


#p.loc

p.year <- ggplot(df, aes(x=plot.pos)) + 
  geom_count(aes(y=sampling.year), color='black', alpha=0.5) + 
  coord_flip() +
  scale_x_discrete(labels = NULL, breaks = df$plot.pos) +
  labs(x = "", y='Year') +
  theme_linedraw() +
  theme(#panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    #panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) +
  ggtitle("Sampling year") #+ theme(legend.position = "none")

#p.year

library(cowplot)
grid::current.viewport()
cowplot::plot_grid(p.dendro, p.blocks, p.genus, p.loc, p.year, 
                   nrow=1, align='h',
                   rel_widths = c(0.7, 2, 0.75, 0.75, 1.2),
                   scale=c(1.08, 1, 1.01, 1.01, 1),
                   labels=c("a", "b", "c", "d", "e"))

### adding gggenes prokka plot

gggenes_df_from_gff_dir <- function(gff_dir) {
  
  gff_files <- list.files(path = gff_dir, full.names = T, pattern = ".gff")
  
  gggenes_df <- data.frame(matrix(nrow = 0 ,ncol = 9))
  
  colnames(gggenes_df) <- c("number","filename_prefix","start","end","gene","strand","direction","type","contig")
  
  for (file in gff_files) {
    print(file)
    
    #set the filename
    filename = tail(unlist(stringr::str_split(file,"/")),1)
    filename_prefix = tail(unlist(stringr::str_split(filename,"[.]gff")),2)[1]
    
    
    #replace any hashes in the name with an underscore and write a temporary file:
    changed_filename_prefix <- gsub("#","_",filename_prefix)
    
    system(glue::glue("sed 's/{filename_prefix}/{changed_filename_prefix}/g' {file} > temp_file.gff"))
    
    #read in the gff file and coerce into a format gggenes likes

    temp_df <- read_delim(file = "temp_file.gff",
                          delim = "\t", escape_double = FALSE, 
                          col_names = FALSE, comment = "#", trim_ws = TRUE) %>%
      rename(contig = X1,method = X2,type = X3, start = X4,end = X5,strand = X6 ,direction = X7,score = X8, details  = X9) %>%
      filter(details != "")  %>% # remove lines where there isn't a 9th column
      mutate(filename_prefix = filename_prefix) %>%
      filter(type == "CDS") %>% # only take CDSs
      # now - if there is ONE  "ID=" string within the details column, create variable "gene" and set as the string immediately after "ID=", but before the next ";". Otherwise set as NA.
      mutate(gene = as.character(purrr::map(details, ~ ifelse(stringr::str_detect(.x, "product=") & length(unlist(stringr::str_split(.x,"product="))) == 2,  unlist(stringr::str_split(unlist(stringr::str_split(.x,"product="))[2],";"))[1] , NA  )))) %>%
      # put in number for gggenes:
      tibble::rownames_to_column(var = "number") %>%
      mutate(number = as.numeric(number)) %>%
      #set the direction and strands:
      mutate(direction = ifelse(direction == "+", 1, -1)) %>%
      mutate(strand = ifelse(direction == 1, "forward","reverse")) %>%
      # now select columns to rbind in:
      select(number,filename_prefix,start,end,gene,strand,direction,type,contig)
    
    
    gggenes_df <- rbind(gggenes_df,temp_df)
    
  }
  
  file.remove("temp_file.gff")
  
  
  return(gggenes_df)
  
  
} # adapted from 'djw533/micro.gen.extra' to replace ID with Name

annot <- gggenes_df_from_gff_dir('/Users/willmatlock/Desktop/gffs')

annot <- annot[!(annot$contig %in% c('NZ_CP073313.1', 'NZ_VAAM01000031.1', 'NZ_UARQ01000013.1')),]
annot <- annot[annot$contig %in% genome_labels,] 


GES.gene.locations <- annot[grepl("GES-5", annot$gene),c("contig","start", "end")]
rownames(GES.gene.locations) <- GES.gene.locations$contig
GES.gene.locations <- GES.gene.locations[-1]

transformAnnot <- function(contig, default_offset=10000){
  offset <- default_offset-GES.gene.locations[contig, c("start")]
  return(offset)
}

#annot$new.start <- apply(annot, MARGIN=1, function(x) as.numeric(x["start"])+transformAnnot(x["contig"]))-10001
#annot$new.end <- apply(annot, MARGIN=1, function(x) as.numeric(x["end"])+transformAnnot(x["contig"]))-10001

annot$contig <- ordered(annot$contig,
                         levels=genome_labels)

library(gggenes)
p.annot <- ggplot(annot, aes(xmin = start, xmax = end, y = contig, fill = gene)) +
  geom_gene_arrow(arrow_body_height = unit(2, "mm"),
                  arrowhead_height = unit(2, "mm"),
                  arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "left") +
  theme_genes()+
  ylab("")+
  xlab("Position (bp)") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(labels=NULL) +
  ggtitle("Prokka annotations")

grid::current.viewport()
cowplot::plot_grid(p.dendro, p.blocks, p.annot,
                   nrow=1, align='h',
                   rel_widths = c(0.7, 2, 2),
                   scale=c(1.08, 1, 1),
                   labels=c("a", "b", "c"))
  

