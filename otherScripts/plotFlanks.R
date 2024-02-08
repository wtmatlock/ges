require(ggplot2)
require(cowplot)
require(gggenes)
require(ggdendro)
require(reshape2)
require(tidyr)
require(vegan)
require(ggdendro)
require(readr)
require(dplyr)
require(stringr)
library(ggfittext)

###
### script parameters
###

pangraph.output <- '/Users/willmatlock/Desktop/thesis-scraps/chapter3/pangraph-output-10000/pangraph-export.gfa.blocks.csv'
gff.directory <- '/Users/willmatlock/Desktop/thesis-scraps/chapter3/gffs' # directory with gff3s for all genomes
ges.block <- 'GWKFGRWJMP' # the block containing blaGES-5
ges.start <- 40 # where blaGES-5 starts in the block
# some contigs might have multiple GES alleles - remove them to avoid repeating the GES block
multiple.alleles <- c('NZ_CP073313.1', 'NZ_VAAM01000031.1', 'NZ_UARQ01000013.1')

###
### generate a gggenes dataframe of contig annotations from NCBI gffs
###

# first, wrangle gffs into gggenes format
gff.files <- list.files(path = gff.directory, full.names = T, pattern = ".gff")
gff.df <- data.frame(matrix(nrow = 0 ,ncol = 9))
colnames(gff.df) <- c("number","filename.prefix","start","end","gene","strand","direction","type","genome")
# some gffs might be empty e.g., not annotated yet
gff.all <- c()
gff.not.empty <- c() 

for (file in gff.files) {
  # set the filename
  filename = tail(unlist(stringr::str_split(file,"/")),1)
  filename.prefix = tail(unlist(stringr::str_split(filename,"[.]gff")),2)[1] 
  gff.all <- c(gff.all, filename.prefix)
  # replace any hashes in the name with an underscore and write a temporary file:
  changed.filename.prefix <- gsub("#","_",filename.prefix)
  system(glue::glue("sed 's/{filename.prefix}/{changed.filename.prefix}/g' {file} > temp_file.gff"))
  #read in the gff file and coerce into a format gggenes likes
  temp.df <- read_delim(file = "temp_file.gff",
                        delim = "\t", escape_double = FALSE, 
                        col_names = FALSE, comment = "#", trim_ws = TRUE) %>%
    rename(contig = X1, method = X2, type = X3, start = X4,end = X5, strand = X6, direction = X7, score = X8, details  = X9) %>%
    filter(details != "") %>% # remove lines where there isn't a 9th column
    mutate(filename.prefix = filename.prefix) %>%
    filter(type %in% c("CDS","gene")) # only take CDSs/genes
  # stop if empty  
  if (nrow(temp.df) > 0){
    temp.df <- temp.df %>%
      group_by(start, end) %>%
      mutate(
        label = if_else(
          (type == "gene" & grepl("gene=([^;]+)", details)),
          stringr::str_extract(details, "(?<=gene=)[^;]+"),
          ifelse(
            (type == "CDS" & grepl("product=([^;]+)", details)),
            stringr::str_extract(details, "(?<=product=)[^;]+"),
            NA
          )
        )
      ) %>%
      # filter results
      # when all labels are NA or none are, take the 'gene' type
      filter(if (all(!is.na(label)) | all(is.na(label))) type=='gene' else all()) %>%
      # otherwise, keep the non-NA gene or CDS label
      filter(if (n()>1) !is.na(label) else all()) %>%
      mutate(label = replace_na(label, "")) %>%
      # put in number for gggenes
      tibble::rownames_to_column(var = "number") %>%
      mutate(number = as.numeric(number)) %>%
      # set the direction and strands
      mutate(direction = ifelse(direction == "+", 1, -1)) %>%
      mutate(strand = as.character(ifelse(direction == 1, "forward","reverse"))) %>%
      # now select columns to rbind in
      select(number, filename.prefix, start, end, label, strand, direction, type, contig)
    gff.not.empty <- c(gff.not.empty, temp.df[1, 'filename.prefix', drop=TRUE])
    gff.df <- rbind(gff.df,temp.df)
  }
}
file.remove("temp_file.gff")
gff.not.empty<- unique(na.omit(gff.not.empty))
gff.empty <- setdiff(gff.all, gff.not.empty)
no.empty <- length(gff.empty)
gff.empty.padding <- data.frame(number = rep(1, no.empty),
                                 filename.prefix = gff.empty,
                                 start = rep(0, no.empty),
                                 end = rep(0, no.empty),
                                 label = rep("GES-5", no.empty),
                                 strand = rep("forward", no.empty),
                                 direction = rep(1, no.empty),
                                 type = rep("gene", no.empty),
                                 contig = gff.empty)
gff.df <- rbind(gff.df, gff.empty.padding)
# check we have all our gffs accounted for
length(unique(gff.df$filename.prefix))==length(gff.all)
# blaGES-5 label
gff.df[grepl('GES-5', gff.df$label),]$label <- "blaGES-5"
# remove contigs with multiple GES alleles, which confuses pangraph
gff.df <- gff.df[!(gff.df$filename.prefix %in% multiple.alleles),]
# add indicator column for plotting later
gff.df$is.integron.finder <- FALSE


###
### add integron finder annotations to the gggenes dataframe
###

integron.finder <- read_delim("/Users/willmatlock/Desktop/thesis-scraps/chapter3/integron-finder.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              comment = "#", trim_ws = TRUE)
integron.finder <- integron.finder[!grepl("ID_integron", integron.finder$`ID_integron`),]
# number annotations for each genome
integron.finder <- integron.finder %>%
  group_by(ID_replicon) %>%
  mutate(number = row_number() + 1000)
integron.finder <- integron.finder[c("number", "ID_replicon", "pos_beg", "pos_end", "annotation", "strand", "strand", "type_elt", "ID_replicon")]
colnames(integron.finder) <- c("number", "filename.prefix", "start", "end", "label", "strand", "direction", "type", "contig")
integron.finder$strand <- ifelse(integron.finder$strand==1, "forward", "backward")
integron.finder$type <- "CDS"
integron.finder <- integron.finder[!grepl("protein", integron.finder$label),]
integron.finder$is.integron.finder <- TRUE
integron.finder$start <- as.numeric(integron.finder$start)
integron.finder$end <- as.numeric(integron.finder$end)
integron.finder$direction <- as.numeric(integron.finder$direction)
gff.df <- rbind(gff.df, integron.finder)

###
### add isfinder annotations to the gggenes dataframe
###

isfinder <- read_delim("/Users/willmatlock/Desktop/thesis-scraps/chapter3/isfinder.tab", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
isfinder <- isfinder[!grepl("#FILE", isfinder$`#FILE`),]
# number annotations for each genome
isfinder <- isfinder %>%
  group_by(SEQUENCE) %>%
  mutate(number = row_number() + 2000)
isfinder$strand <- ifelse(isfinder$STRAND=="+", "forward", "reverse")
isfinder$STRAND <- ifelse(isfinder$strand=="forward", 1, -1)
isfinder <- isfinder[c("number", "SEQUENCE", "START", "END", "GENE", "strand", "STRAND", "DATABASE", "SEQUENCE")]
colnames(isfinder) <- c("number", "filename.prefix", "start", "end", "label", "strand", "direction", "type", "contig")
isfinder$is.integron.finder <- FALSE
isfinder$is.isfinder <- TRUE
gff.df$is.isfinder <- FALSE
isfinder$direction <- as.numeric(isfinder$direction)
isfinder$start <- as.numeric(isfinder$start)
isfinder$end <- as.numeric(isfinder$end)
gff.df <- rbind(gff.df, isfinder)

# filter out NCBI transposase annotations in favour of ISfinder

gff.df <- gff.df[!(gff.df$is.isfinder==FALSE & grepl("transposase", gff.df$label)),]


###
### add pangraph blocks to the gggenes dataframe
###

# genome blocks approach
genome.blocks <- read.csv(pangraph.output, header=T, stringsAsFactors = F)
genome.blocks <- genome.blocks[!(genome.blocks$genome %in% multiple.alleles),]
genome.blocks$forward <- ifelse(genome.blocks$strand=="+", TRUE, FALSE)
block.counts <- table(genome.blocks$block)
blocks.which.need.colours <- names(block.counts)[which(block.counts>1)]
genome.blocks$block.coloured <- sapply(genome.blocks$block,
                                       function(x) ifelse(x %in% blocks.which.need.colours,
                                                          x,
                                                          "_other"))
block.colours <- unique(genome.blocks$colour)
names(block.colours) <- unique(genome.blocks$block.coloured)
ges.block.locations <- genome.blocks[which(genome.blocks$block==ges.block),c("start", "end")]
rownames(ges.block.locations) <- genome.blocks[which(genome.blocks$block==ges.block),"genome"]
transformAnnotationsBlocks <- function(genome, default_offset=10000){
  offset <- default_offset - ges.block.locations[genome, c("start")]
  return(offset)
}

genome.blocks$new.start <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["start"])+transformAnnotationsBlocks(x["genome"]))-10001-ges.start
genome.blocks$new.end <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["end"])+transformAnnotationsBlocks(x["genome"]))-10001 - ges.start
# check agreement between datasets, modulo contigs removed deliberately
setdiff(genome.blocks$genome, gff.df$filename.prefix)
setdiff(gff.df$filename.prefix, genome.blocks$genome)
# get all the genome paths (in terms of blocks)
genome.paths <- sapply(unique(genome.blocks$genome), function(x) paste(genome.blocks[which(genome.blocks$genome==x), "block"], collapse=","))
genome.paths <- sort(genome.paths, decreasing = TRUE)
# subset to unique paths, keeping one representative genome for each
genome.path.reps <- names(genome.paths)[!duplicated(genome.paths)]
genome.blocks.unique <- genome.blocks[which(genome.blocks$genome %in% genome.path.reps),]
# store the number of examples of each
genome.blocks.unique$genome.path <- genome.paths[genome.blocks.unique$genome]
genome.blocks.unique$n.reps <- sapply(genome.blocks.unique$genome.path,
                                      function(x) table(genome.paths)[x])# all genomes for each rep
genome.blocks.unique$all.reps <- sapply(genome.blocks.unique$genome.path,
                                        function(x) names(which(genome.paths==x)))
genome.blocks.unique$genome.path.name <- paste0("Type", as.numeric(as.factor(genome.blocks.unique$genome.path )))
genome.blocks.unique$genome.n <- sapply(genome.blocks.unique$n.reps, function(x) ifelse(x==1,
                                                                                        "", paste0("n=", x)))

###
### flipping and trimming the gffs
###

gff.df.trim <- gff.df %>%
  filter(filename.prefix %in% genome.path.reps) %>%
  group_by(contig) %>%
  # first, we shift genes so blaGES-5 starts at 0
  mutate(gene.length = end - start,
         ges.length = gene.length[label=='blaGES-5'],
         ges.shift = start[label=='blaGES-5'],
         start = ifelse(label=='blaGES-5', 0, start - ges.shift),
         end = start + gene.length,
         # then, if blaEGS-5 is reversed, we flip everything around
         is.ges.reversed = ifelse(direction[label=='blaGES-5']==-1, TRUE, FALSE),
         new.start = ifelse(is.ges.reversed==TRUE, -end + ges.length, start),
         new.end= ifelse(is.ges.reversed==TRUE, -start + ges.length, end),
         new.direction = ifelse(is.ges.reversed, -as.numeric(direction), as.numeric(direction))) %>%
  # filter for genes within the 10kbp flanking region
  filter(!(new.start < -10000)) %>%
  filter(!(new.end > 10000))

###
### clustering contigs by gene cassette blocks
###

genome.blocks.unique %>% 
  group_by(genome) %>%
  summarise(reps = lengths(all.reps)) %>%
  summary()
  
m <- acast(genome ~ block, data=genome.blocks.unique, fill=0, fun.aggregate=length)
# count blocks downstream of GES-5
block.counts.df <- as.data.frame(block.counts)
colnames(block.counts.df) <- c("block", "count")
# get integron blocks
cassette.boundaries <- gff.df.trim %>%
  filter(is.integron.finder==TRUE) %>%
  filter(label %in% c("attC", "attI_1", "attI_3")) %>%
  group_by(contig) %>%
  mutate(n.att = n()) %>%
  filter(any(n.att >= 2)) %>%
  group_by(contig) %>%
  summarise(att.start = min(new.start),
            att.end = max(new.end))

block.repertoire <- merge(genome.blocks.unique, cassette.boundaries, by.x='genome', by.y='contig', all.y=TRUE)
block.repertoire <- block.repertoire %>%
  group_by(genome) %>%
  filter(new.start >= att.start) %>%
  filter(new.end <= att.end)
block.repertoire <- unique(block.repertoire$block)
m.filtered <- m[,block.repertoire]
m.filtered <- m.filtered[rowSums(m.filtered[])>0,]
m.dist <- vegdist(m.filtered, method="jaccard") # jaccard distances based on block presence/absence
# clustering
dendro <- hclust(m.dist, method="complete")
dendro_order <- order.dendrogram(as.dendrogram(dendro))
genome_labels <- dendro_data(dendro)$labels$label
genome.blocks.unique$genome.ordered <- ordered(genome.blocks.unique$genome,
                                               levels=genome_labels)


# now order gffs by the dendrogram
gff.df.trim$genome.ordered <- ordered(gff.df.trim$contig,
                                               levels=genome_labels)

# separate out promoter annotations for plotting later
promoters <- gff.df.trim[grepl("Pc|Pint", gff.df.trim$label),][c("contig", "label")]
gff.df.trim <- gff.df.trim[!grepl("Pc|Pint", gff.df.trim$label),]

# position of representative genomes from bottom to top for plotting metadata later
position <- cbind(genome_labels, length(genome_labels):1)
plot.placement <- merge(genome.blocks.unique, position, by.x='genome.ordered', by.y='genome_labels', all.x=TRUE)
plot.placement <- plot.placement[,c("genome.ordered", "all.reps", "V2")]
colnames(plot.placement) <- c("genome.ordered", "all.reps", "pos")
library(tidyr)
plot.placement <- unnest(plot.placement, cols=all.reps)
plot.placement <- unique(plot.placement)
plot.placement$pos <- as.numeric(plot.placement$pos)

### block distribution


blocks.dist <- genome.blocks %>%
  group_by(block) %>%
  mutate(n.block = n(),
         length = abs(start-end)) %>%
  select(block, n.block, length) %>%
  distinct()

summary(blocks.dist$length)
sum(blocks.dist$length)
sum(blocks.dist[blocks.dist$n.block==1,]$length)

###
### plot pangraph blocks with integron finder and gff annotations
###

# lighten block colours
library(colorspace)
new.block.colours <- desaturate(lighten(block.colours, 0.7), 0.3)

gff.df.trim <- gff.df.trim[!is.na(gff.df.trim$genome.ordered),]

# add integrase colours
library(readr)
int.clusters <- read_csv("/Users/willmatlock/Desktop/thesis-scraps/chapter3/integrase-snps/integrase-clusters.csv")
colnames(int.clusters) <- c("contig", "cluster")

gff.df.trim <- merge(gff.df.trim, int.clusters, by.x='filename.prefix', by.y="contig", all.x=TRUE)
gff.df.trim[!grepl("intI", gff.df.trim$label),]$cluster <- NA
gff.df.trim[gff.df.trim$is.integron.finder==FALSE,]$cluster <- NA

p.flanks <- ggplot() +
 geom_gene_arrow(data = genome.blocks.unique,
                  aes(xmin = new.start, xmax = new.end, forward = forward, y = genome.ordered, fill = block.coloured),
                  colour = NA,
                  arrow_body_height = unit(3, "mm"),
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(0, "mm")) +
  theme_genes() +
  scale_fill_manual(values=new.block.colours) +
  ylab("") +
  theme(legend.position = "none") +
  scale_y_discrete(breaks=genome.blocks.unique$genome.ordered, labels=genome.blocks.unique$genome.n) +
  geom_gene_arrow(data=gff.df.trim[gff.df.trim$is.integron.finder==FALSE,], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="black", fill=NA,
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim[gff.df.trim$is.integron.finder==TRUE & grepl("attC", gff.df.trim$label),], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="#e41a1c", fill='#e41a1c',
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim[gff.df.trim$is.integron.finder==TRUE & grepl("attI", gff.df.trim$label),], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="#ff7f00", fill='#ff7f00',
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim[gff.df.trim$is.integron.finder==TRUE & grepl("int", gff.df.trim$label) &gff.df.trim$cluster==1,], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="black", fill='#377eb8',
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim[gff.df.trim$is.integron.finder==TRUE & grepl("int", gff.df.trim$label) &gff.df.trim$cluster==2,], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="black", fill='#4daf4a',
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm"))

###
### plot integrase promoters heatmap
###

promoters <- merge(promoters, plot.placement, by.x='contig', by.y='all.reps', all.y=TRUE)


promoters$label <- factor(promoters$label, levels = c("Pint_1", "Pc_1"))
promoters$value <- 1
promoters <- promoters %>%
  pivot_wider(names_from = label, values_from = value, values_fill = 0) %>%
  mutate(Pint_1 = ifelse(Pint_1, "1", "0"), Pc_1 = ifelse(Pc_1, "1", "0"))
promoters <- promoters[c("contig", "pos", "Pint_1", "Pc_1")]
promoters <- melt(promoters, id=c("contig", "pos"))
promoters$pos <- as.numeric(promoters$pos)
promoters$plot.pos <- abs(promoters$pos - max(promoters$pos)) +1

p.promoters <- ggplot(promoters, aes(x = variable, y = plot.pos, fill=value)) +
  geom_tile(colour="black", size=0.25) +
  scale_y_discrete(labels = NULL, breaks = promoters$plot.pos) +
  theme(legend.position = "none") +
  labs(x = "", y='') +
  scale_x_discrete(labels=c("Pint", "Pc")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) +
  scale_fill_manual(values=c("white", "black"))
  
###
### plot intI1 SNVs
###

int.snv <- read_table("/Users/willmatlock/Desktop/thesis-scraps/chapter3/integrase-snps/C1/C1.snv.fasta.tab", 
                           col_names = FALSE)
colnames(int.snv) <- c("contig", "snv-profile")
int.snv$contig <- gsub("-integrase", "", int.snv$contig)
int.snv <- merge(int.snv, plot.placement, by.x='contig', by.y='all.reps', all.y=TRUE)
int.snv$pos <- as.numeric(int.snv$pos)
int.snv$plot.pos <- abs(int.snv$pos - max(int.snv$pos)) +1
int.snv[int.snv$contig %in% c("NZ_CP081345.1", "NZ_CP107042.1"),]$`snv-profile` <- "CTGA/CTGC"

p.int <- ggplot(int.snv, aes(x = "", y = plot.pos, fill=`snv-profile`)) +
  geom_tile(colour="black", size=0.25) +
  scale_x_discrete(labels=c("intI1")) +
  scale_y_discrete(labels = NULL, breaks = int.snv$plot.pos) +
  #theme(legend.position = "none") +
  labs(x = "", y='') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) +
  scale_fill_brewer(palette="Pastel2")

## plot ges SNVs

ges.snv <- read_table("/Users/willmatlock/Desktop/thesis-scraps/chapter3/ges-snps/ges.snps.tab", 
                      col_names = FALSE)
colnames(ges.snv) <- c("contig", "snv-profile")
ges.snv <- merge(ges.snv, plot.placement, by.x='contig', by.y='all.reps', all.y=TRUE)
ges.snv$pos <- as.numeric(ges.snv$pos)
ges.snv$plot.pos <- abs(ges.snv$pos - max(ges.snv$pos)) +1

p.ges <- ggplot(ges.snv, aes(x = "", y = plot.pos, fill=`snv-profile`)) +
  geom_tile(colour="black", size=0.25) +
  #geom_text(aes(label=`snv-profile`), size=2.5) +
  scale_x_discrete(labels=c("blaGES-5")) +
  scale_y_discrete(labels = NULL, breaks = ges.snv$plot.pos) +
  #theme(legend.position = "none") +
  labs(x = "", y='') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) +
  scale_fill_brewer(palette="Pastel1")

# plot contig origin - chromosome, plasmids, NA

mlst <- read_delim("/Users/willmatlock/Desktop/thesis-scraps/chapter3/mlst-summary.tab", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE, col_names = FALSE)
mlst$contig <- word(mlst$X1, 3, sep="/")
mlst$contig <- gsub(".fasta", "", mlst$contig)
mlst <- mlst[mlst$contig %in% unique(genome.blocks$genome),]
mlst$location <- ifelse(mlst$X2!="-", "chromosome", NA)
mlst$location.info <- paste0(mlst$X2, mlst$X3)

loc <- mlst[c("contig", "location", "location.info")]

plasmidfinder <- read_delim("/Users/willmatlock/Desktop/thesis-scraps/chapter3/plasmidfinder_summary.tab", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
plasmidfinder$contig <- word(plasmidfinder$`#FILE`, 3, sep="/")
plasmidfinder <- plasmidfinder[plasmidfinder$contig %in% unique(genome.blocks$genome),]
plasmidfinder[3:7] <- data.frame(lapply(plasmidfinder[3:7], function(x) { gsub("\\.", "", x) }))
plasmidfinder[3:7] <- Map(function(n, x) replace(x, x != "", n), names(plasmidfinder[3:7]), plasmidfinder[3:7])
plasmidfinder$location.info <- do.call(paste, c(plasmidfinder[3:7], sep=""))
plasmidfinder$location <- ifelse(plasmidfinder$location.info!="", "plasmid", NA)
plasmidfinder <- plasmidfinder[c("contig", "location", "location.info")]

loc <- rbind(loc, plasmidfinder)
loc <- loc[!is.na(loc$location),]

### additional mobtype replicons

mobtyper <- read_table("/Users/willmatlock/Desktop/thesis-scraps/chapter3/mobtyper-summary.tab")
mobtyper <- mobtyper[!grepl("sample_id", mobtyper$sample_id),]
mobtyper <- mobtyper[!grepl("sample_id", mobtyper$sample_id),]
mobtyper <- mobtyper[c("sample_id", "rep_type(s)", "predicted_mobility")]

loc <- merge(loc, mobtyper, by.x="contig", by.y="sample_id", all.y=TRUE)
loc[is.na(loc$location)&!is.na(loc$`rep_type(s)`)&loc$`rep_type(s)`!="-",]$location <- "plasmid"
loc <- loc[!(loc$contig %in% multiple.alleles),]

summary(as.factor(loc[loc$location=="plasmid",]$predicted_mobility))

# import metadata

library(readr)
library(dplyr)
library(stringr)

metadata <- read_delim("/Users/willmatlock/Desktop/thesis-scraps/chapter3/ges/data/NCBI-GES.tsv", 
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

who_region <- read_csv("/Users/willmatlock/Desktop/thesis-scraps/chapter3/ges/data/who-region.csv")
who_region <- who_region[,c("name", 'region')]
df <- merge(df, who_region, by.x='country', by.y='name', all.x=TRUE)
df$region <- ifelse(df$country=='Czech Republic', "Europe",
                    ifelse(df$country=="Russia", "Europe",
                           ifelse(df$country=="Taiwan", "Asia",
                                  ifelse(df$country=="United Kingdom", 'Europe',
                                         ifelse(df$country=="USA", "Americas", df$region))))) 

df <- merge(df, loc, by.x="genome.ordered", by.y="contig", all.x=TRUE)


p.location <- ggplot(df, aes(x = "", y = plot.pos, fill=location)) +
  geom_tile(colour="black", size=0.25) +
  scale_x_discrete(labels=c("Replicon")) +
  scale_y_discrete(labels = NULL, breaks = df$plot.pos) +
  #theme(legend.position = "none") +
  labs(x = "", y='') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) +
  scale_fill_brewer(palette="Pastel1")
  
p.genus <- ggplot(df, aes(x=plot.pos)) + 
  geom_bar(aes(y=..count.., fill=genus), position='fill',size=0.25, color='black', width = 1) + 
  coord_flip() +
  scale_x_discrete(labels = NULL, breaks = df$plot.pos) +
  scale_y_continuous(breaks = c(0.0, 1.0)) +
  labs(x = "", y='') +
  scale_fill_brewer("", palette = "Set3") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) + 
  #theme(legend.position = "none")

p.loc <- ggplot(df, aes(x=plot.pos)) + 
  geom_bar(aes(y=..count.., fill=region), position='fill',size=0.25, color='black', width = 1) + 
  coord_flip() +
  scale_x_discrete(labels = NULL, breaks = df$plot.pos) +
  scale_y_continuous(breaks = NULL) +
  labs(x = "", y='') +
  scale_fill_brewer("", palette="Set2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) + 
  #theme(legend.position = "none")

p.year <- ggplot(df, aes(x=plot.pos)) + 
  geom_count(aes(y=sampling.year), color='black', alpha=0.5) + 
  coord_flip() +
  scale_x_discrete(labels = NULL, breaks = df$plot.pos) +
  labs(x = "", y='') +
  theme_linedraw() +
  theme(#panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    #panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) + 
  theme(legend.position = "none")


library(cowplot)
grid::current.viewport()
cowplot::plot_grid(p.flanks,
                   p.promoters,
                   p.int,
                   p.ges,
                   p.location,
                   p.genus,
                   p.loc,
                   p.year, 
                   nrow=1, align='h',
                   rel_widths = c(7, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 3),
                   scale=c(1, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1),
                   labels = c("a", "b", "c", "d", "e", "f", "g", "h"))

###
# subset plot to integrase and gene cassettes 
###

gff.df.trim.trim <- gff.df.trim[gff.df.trim$new.start>=-3000 & gff.df.trim$new.end<=5000,]

### manually curating gene labels

gff.df.trim.trim$new.label <- ""
gff.df.trim.trim[grepl("GES-5", gff.df.trim.trim$label),]$new.label <- "blaGES-5"
gff.df.trim.trim[grepl("intI", gff.df.trim.trim$label),]$new.label <- ""
gff.df.trim.trim[grepl("AAC\\(6'\\)-Ib4", gff.df.trim.trim$label),]$new.label <- "aac(6')-Ib4"
gff.df.trim.trim[grepl("DUF1010", gff.df.trim.trim$label),]$new.label <- "DUF1010"
gff.df.trim.trim[grepl("DUF3788", gff.df.trim.trim$label),]$new.label <- "DUF3788"
#gff.df.trim.trim[grepl("IS110", gff.df.trim.trim$label),]$new.label <- "IS110"
gff.df.trim.trim[grepl("AAC\\(6'\\)-Ib", gff.df.trim.trim$label),]$new.label <- "aac(6')-Ib"
gff.df.trim.trim[grepl("QacE", gff.df.trim.trim$label),]$new.label <- "qacE"
gff.df.trim.trim[grepl("sul1", gff.df.trim.trim$label),]$new.label <- "sul1"
gff.df.trim.trim[grepl("aph\\(3'\\)-XV", gff.df.trim.trim$label),]$new.label <- "aph(3')-XV"
gff.df.trim.trim[grepl("AAC\\(6'\\)-IIa", gff.df.trim.trim$label),]$new.label <- "aac(6')-IIa"
gff.df.trim.trim[grepl("AadA", gff.df.trim.trim$label),]$new.label <- "aadA"
gff.df.trim.trim[grepl("aac\\(6'\\)-31", gff.df.trim.trim$label),]$new.label <- "aac(6')-31"
gff.df.trim.trim[grepl("ANT\\(3''\\)-Ia", gff.df.trim.trim$label),]$new.label <- "aadA4"
gff.df.trim.trim[grepl("qacG2", gff.df.trim.trim$label),]$new.label <- "qacG2"
gff.df.trim.trim[grepl("QacF", gff.df.trim.trim$label),]$new.label <- "qacF"
gff.df.trim.trim[grepl("dfrA5", gff.df.trim.trim$label),]$new.label <- "dfrA5"
gff.df.trim.trim[grepl("ere\\(A\\)", gff.df.trim.trim$label),]$new.label <- "ere(A)"
gff.df.trim.trim[grepl("aadB", gff.df.trim.trim$label),]$new.label <- "aadB"
gff.df.trim.trim[grepl("aac\\(6'\\)", gff.df.trim.trim$label),]$new.label <- "aac(6')"
gff.df.trim.trim[grepl("AAC\\(6'\\)-Ib3", gff.df.trim.trim$label),]$new.label <- "aac(6')-Ib3"
gff.df.trim.trim[grepl("dfrA7", gff.df.trim.trim$label),]$new.label <- "dfrA7"
gff.df.trim.trim[grepl("hypothetical protein", gff.df.trim.trim$label),]$new.label <- "hyp."
gff.df.trim.trim[grepl("aadA1", gff.df.trim.trim$label),]$new.label <- "aadA1"
gff.df.trim.trim[grepl("OXA-2", gff.df.trim.trim$label),]$new.label <- "blaOXA-2"
gff.df.trim.trim[grepl("dfrA21", gff.df.trim.trim$label),]$new.label <- "dfrA21"
gff.df.trim.trim[grepl("OXA-1", gff.df.trim.trim$label),]$new.label <- "blaOXA-1"
gff.df.trim.trim[grepl("OXA-10", gff.df.trim.trim$label),]$new.label <- "blaOXA-17"
gff.df.trim.trim[grepl("qac", gff.df.trim.trim$label),]$new.label <- "qac"
gff.df.trim.trim[grepl("AAC\\(6'\\)-Ia", gff.df.trim.trim$label),]$new.label <- "aac(6')-Ia"
gff.df.trim.trim[grepl("DM13", gff.df.trim.trim$label),]$new.label <- "DM13"
gff.df.trim.trim[grepl("QnrVC4", gff.df.trim.trim$label),]$new.label <- "QnrVC4"
gff.df.trim.trim[grepl("APH\\(3'\\)-XV", gff.df.trim.trim$label),]$new.label <- "aph(3')-XV"
#gff.df.trim.trim[grepl("IS6100", gff.df.trim.trim$label),]$new.label <- "IS6100"
gff.df.trim.trim[grepl("RepB", gff.df.trim.trim$label),]$new.label <- "repB"
gff.df.trim.trim[grepl("MobB", gff.df.trim.trim$label),]$new.label <- "mobB"
gff.df.trim.trim[grepl("mobilization protein", gff.df.trim.trim$label),]$new.label <- "mob. protein"
gff.df.trim.trim[grepl("MobL", gff.df.trim.trim$label),]$new.label <- "mobA/mobL"
gff.df.trim.trim[gff.df.trim.trim$is.isfinder==TRUE,]$new.label <-gff.df.trim.trim[gff.df.trim.trim$is.isfinder==TRUE,]$label 

p.flanks.subset <- ggplot() +
  geom_gene_arrow(data = genome.blocks.unique[genome.blocks.unique$new.start>=-3000 & genome.blocks.unique$new.end<=5000,],
                  aes(xmin = new.start, xmax = new.end, forward = forward, y = genome.ordered, fill = block.coloured),
                  colour = NA,
                  arrow_body_height = unit(3, "mm"),
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(0, "mm")) +
  theme_genes() +
  scale_fill_manual(values=new.block.colours) +
  ylab("") +
  theme(legend.position = "none") +
  scale_y_discrete(breaks=genome.blocks.unique$genome.ordered, labels=genome.blocks.unique$genome.n) +
  geom_gene_arrow(data=gff.df.trim.trim[gff.df.trim.trim$is.integron.finder==FALSE,], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="black", fill=NA,
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim.trim[gff.df.trim.trim$is.integron.finder==TRUE & grepl("attC", gff.df.trim.trim$label),], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="#e41a1c", fill='#e41a1c',
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim.trim[gff.df.trim.trim$is.integron.finder==TRUE & grepl("attI", gff.df.trim.trim$label),], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="#ff7f00", fill='#ff7f00',
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim.trim[gff.df.trim.trim$is.integron.finder==TRUE & grepl("int", gff.df.trim.trim$label) &gff.df.trim.trim$cluster==1,], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="black", fill='#377eb8',
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim.trim[gff.df.trim.trim$is.integron.finder==TRUE & grepl("int", gff.df.trim.trim$label) &gff.df.trim.trim$cluster==2,], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="black", fill='#4daf4a',
                  arrow_body_height = unit(1.8, "mm"),
                  arrowhead_height = unit(1.8, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_label(data=gff.df.trim.trim[gff.df.trim.trim$is.integron.finder==FALSE,], 
                  aes(xmin = new.start, xmax = new.end, label=new.label, y = genome.ordered), 
                  align = "left", reflow=TRUE, grow=TRUE, height = grid::unit(1, "mm")) 

### export df as table S1

table.s1 <- df
table.s1 <- table.s1[c("genome.ordered", "pos", "location", "location.info", "rep_type(s)", "predicted_mobility", "genus", "region",  "sampling.year")]
colnames(table.s1) <- c("contig", "plot.position", "replicon", "plasmid.finder.or.mlst", "mob.typer", "predicted_mobility", "ncbi.genus", "iso.region", "sampling.year")
write_csv(table.s1, '~/Desktop/tableS1.csv')

###
### subset plot to intI1 and SNV profiles
###
  
intI1.genomes <- int.snv %>%
  filter(!is.na(`snv-profile`)) %>%
  group_by(`snv-profile`) %>%
  mutate(freq=n()) %>%
  arrange(desc(freq))

genome.blocks.unique <- genome.blocks.unique[genome.blocks.unique$genome %in% intI1.genomes$contig,]
genome.blocks.unique$genome.ordered <- ordered(genome.blocks.unique$genome,
                                               levels=intI1.genomes$contig)
gff.df.trim <- gff.df.trim[gff.df.trim$contig %in% intI1.genomes$contig,]
gff.df.trim$genome.ordered <- ordered(gff.df.trim$contig,
                                               levels=intI1.genomes$contig)

p.flanks.intI1 <- ggplot() +
  geom_gene_arrow(data = genome.blocks.unique,
                  aes(xmin = new.start, xmax = new.end, forward = forward, y = genome.ordered, fill = block.coloured),
                  colour = NA,
                  arrow_body_height = unit(6, "mm"),
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(0, "mm")) +
  theme_genes() +
  scale_fill_manual(values=new.block.colours) +
  ylab("") +
  theme(legend.position = "none") +
  scale_y_discrete(breaks=genome.blocks.unique$genome.ordered, labels=genome.blocks.unique$genome.n) +
  geom_gene_arrow(data=gff.df.trim[gff.df.trim$is.integron.finder==FALSE,], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="black", fill=NA,
                  arrow_body_height = unit(3, "mm"),
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim[gff.df.trim$is.integron.finder==TRUE & grepl("attC", gff.df.trim$label),], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="#e41a1c", fill='#e41a1c',
                  arrow_body_height = unit(3, "mm"),
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim[gff.df.trim$is.integron.finder==TRUE & grepl("attI", gff.df.trim$label),], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="#ff7f00", fill='#ff7f00',
                  arrow_body_height = unit(3, "mm"),
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(0.5, "mm")) +
  geom_gene_arrow(data=gff.df.trim[gff.df.trim$is.integron.finder==TRUE & grepl("int", gff.df.trim$label) &gff.df.trim$cluster==1,], 
                  aes(xmin = new.start, xmax = new.end, forward = new.direction, y = genome.ordered), 
                  size=0.4, colour="black", fill='#377eb8',
                  arrow_body_height = unit(3, "mm"),
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(0.5, "mm"))

p.int.intI1 <- ggplot(intI1.genomes[intI1.genomes$contig==intI1.genomes$genome.ordered,], aes(x = "", y=(1:25), fill=`snv-profile`)) +
  geom_tile(colour="black", size=0.25) +
  geom_text(aes(label=`snv-profile`), size=2.5) +
  scale_x_discrete(labels=c("intI1")) +
  scale_y_discrete(labels = NULL, breaks = int.snv$plot.pos) +
  theme(legend.position = "none") +
  labs(x = "", y='') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) +
  scale_fill_brewer(palette="Pastel2")


library(cowplot)
grid::current.viewport()
cowplot::plot_grid(p.flanks.intI1,
                   p.int.intI1,
                   nrow=1, align='h',
                   rel_widths = c(7, 1),
                   scale=c(1, 1.03),
                   labels = c("a", "b"))

###
### gene cassette distribution
###

# of the 104 genomes, 102 have gff annotations 
# the 104 genomes represent 80 unique block paths

reps <- genome.blocks.unique %>%
  mutate(reps = lengths(all.reps)) %>%
  select(genome, reps) %>%
  distinct()

df.cassette <- gff.df.trim %>%
  arrange(contig, new.start) %>%
  group_by(contig) %>%
  mutate(prev_annotation = lag(label),
         next1_annotation = lead(label),
         next2_annotation = lead(label, n=2),
         next3_annotation = lead(label, n=3),
         next4_annotation = lead(label, n=4),
         next5_annotation = lead(label, n=5)) %>%
  mutate(gene.cassette = ifelse(prev_annotation %in% c("attC", "attI_1") & next1_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                label, ifelse(
                                  prev_annotation %in% c("attC", "attI_1") & next2_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                  paste(label, lead(label), sep="|" ), ifelse(
                                    prev_annotation %in% c("attC", "attI_1") & next3_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                    paste(paste(label, lead(label), sep="|" ), lead(label, n=2), sep="|"), ifelse(
                                      prev_annotation %in% c("attC", "attI_1") & next4_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                      paste(paste(paste(label, lead(label), sep="|" ), lead(label, n=2), sep="|"), lead(label, n=3), sep="|"), NA))))) %>%
  ungroup() %>%
  select(contig, gene.cassette) %>%
  filter(!is.na(gene.cassette)) %>%
  mutate(n.genes = 1 + str_count(gene.cassette,"\\|"))

df.cassette <- merge(df.cassette, reps, by.x='contig', by.y='genome', all.x=TRUE)

df.cassette <- df.cassette %>%
  group_by(gene.cassette) %>%
  mutate(total = sum(reps)) %>%
  select(gene.cassette, total, n.genes) %>%
  distinct()

df.cassette$gene.cassette <- paste0(paste0("|", df.cassette$gene.cassette), "|")

summary(df.cassette$n.genes)

# write.csv(df.cassette, '~/Desktop/cassette_distribution.csv')

###
### gene cassette positions
###

# subset df.cassette to only contigs with integrase and two attachment sites

reps <- genome.blocks.unique %>%
  mutate(reps = lengths(all.reps)) %>%
  select(genome, reps) %>%
  distinct()

df.position <- gff.df.trim %>%
  arrange(contig, new.start) %>%
  group_by(contig) %>%
  filter(any(is.integron.finder==TRUE & (label=="intI" | label %in% c("attI_1", "attI_3")))) %>%
  mutate(prev_annotation = lag(label),
         next1_annotation = lead(label),
         next2_annotation = lead(label, n=2),
         next3_annotation = lead(label, n=3),
         next4_annotation = lead(label, n=4),
         next5_annotation = lead(label, n=5)) %>%
  mutate(gene.cassette = ifelse(prev_annotation %in% c("attC", "attI_1") & next1_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                label, ifelse(
                                  prev_annotation %in% c("attC", "attI_1") & next2_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                  paste(label, lead(label), sep="|" ), ifelse(
                                    prev_annotation %in% c("attC", "attI_1") & next3_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                    paste(paste(label, lead(label), sep="|" ), lead(label, n=2), sep="|"), ifelse(
                                      prev_annotation %in% c("attC", "attI_1") & next4_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                      paste(paste(paste(label, lead(label), sep="|" ), lead(label, n=2), sep="|"), lead(label, n=3), sep="|"), NA))))) %>%
  filter(!is.na(gene.cassette)) %>%
  mutate(order = row_number()) %>%
  select(contig, gene.cassette, order)

df.position <- merge(df.position, reps, by.x='contig', by.y='genome', all.x=TRUE)

df.position <- df.position %>%
  group_by(contig, gene.cassette, order) %>%
  slice(rep(row_number(), reps)) %>%
  ungroup() %>%
  select(gene.cassette, order) %>%
  group_by(gene.cassette) %>%
  summarise(median = median(order),
         min = min(order),
         max = max(order))

df.position$gene.cassette <- paste0(paste0("|", df.position$gene.cassette), "|")

df.cassette <- merge(df.cassette, df.position, by.x='gene.cassette', by.y='gene.cassette', all.x=TRUE)
df.cassette <- df.cassette[order(-df.cassette$total, df.cassette$gene.cassette),]

write.csv(df.cassette, '~/Desktop/cassette_distribution.csv')

###
### gene array distribution
###

df.array <- gff.df.trim %>%
  arrange(contig, new.start) %>%
  group_by(contig) %>%
  filter(any(is.integron.finder==TRUE & (label=="intI" | label %in% c("attI_1", "attI_3")))) %>%
  mutate(prev_annotation = lag(label),
         next1_annotation = lead(label),
         next2_annotation = lead(label, n=2),
         next3_annotation = lead(label, n=3),
         next4_annotation = lead(label, n=4),
         next5_annotation = lead(label, n=5)) %>%
  mutate(gene.cassette = ifelse(prev_annotation %in% c("attC", "attI_1") & next1_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                label, ifelse(
                                  prev_annotation %in% c("attC", "attI_1") & next2_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                  paste(label, lead(label), sep="|" ), ifelse(
                                    prev_annotation %in% c("attC", "attI_1") & next3_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                    paste(paste(label, lead(label), sep="|" ), lead(label, n=2), sep="|"), ifelse(
                                      prev_annotation %in% c("attC", "attI_1") & next4_annotation %in% c("attC") & !(label %in% c("attC", "attI_1")),
                                      paste(paste(paste(label, lead(label), sep="|" ), lead(label, n=2), sep="|"), lead(label, n=3), sep="|"), 
                                      ifelse(is.integron.finder==TRUE, label, NA)))))) %>%
  filter(!is.na(gene.cassette)) %>%
  mutate(order = row_number()) %>%
  select(contig, gene.cassette, order, cluster) %>%
  group_by(contig) %>%
  mutate(array = paste(gene.cassette, collapse="|")) %>%
  select(contig, cluster, array) %>%
  filter(!is.na(cluster)) %>%
  group_by(array) %>%
  mutate(n=n())


write.csv(df.array, '~/Desktop/array_distribution.csv')

