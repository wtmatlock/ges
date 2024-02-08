require(ggplot2)
require(cowplot)
require(gggenes)
require(ggdendro)
require(reshape2)
require(tidyr)
require(vegan) # for Jaccard distance
require(ggdendro)
require(readr)
require(dplyr)
require(stringr)

### script inputs

pangraph.output <- '/Users/willmatlock/Desktop/pangraph-output-10000/pangraph-export.gfa.blocks.csv'
gff.directory <- '/Users/willmatlock/Desktop/gffs' # directory with gff3s for all genomes
ges.block <- 'GWKFGRWJMP' # the block containing blaGES-5
ges.start <- 40 # where blaGES-5 starts in the block
# some contigs might have multiple GES alleles, likely due to poor assembly - remove them to avoid repeating the GES block
multiple.alleles <- c('NZ_CP073313.1', 'NZ_VAAM01000031.1', 'NZ_UARQ01000013.1')

### first, wrangle gffs into gggenes format

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

# general formatting of CDS labels - will probably need manual curation!
gff.df[gff.df$type=="CDS"&gff.df$label!="blaGES-5",]$label <- str_extract_all(gff.df[gff.df$type=="CDS"&gff.df$label!="blaGES-5",]$label, "\\b[A-Z][a-zA-Z0-9_-]{2,}\\b") %>%
  lapply(function(x) ifelse(length(x) > 0, tail(x, n = 1), "")) %>%
  unlist()

# remove contigs with multiple GES alleles, which confuses pangraph
gff.df <- gff.df[!(gff.df$filename.prefix %in% multiple.alleles),]

# genome blocks approach
genome.blocks <- read.csv(pangraph.output, header=T, stringsAsFactors = F)
genome.blocks <- genome.blocks[!(genome.blocks$genome %in% multiple.alleles),]
genomes <- unique(genome.blocks$genome)

plot.df <- gff.df[gff.df$filename.prefix %in% genomes,]

lengths <- read_delim("/Users/willmatlock/Desktop/ges/data/lengths.tsv", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
colnames(lengths) <- c("contig", "length")
lengths$contig <- stringr::word(lengths$contig, 1)
lengths <- lengths[lengths$contig %in% genomes,]

df <- read_delim("~/Desktop/ges/data/NCBI-GES.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
df <- df[df$Contig %in% genomes,]
df$genus <- word(df$`#Scientific name`, 1)
df$country <- word(df$Location, 1, sep=":")
df$sampling.year <- substr(df$`Collection date`, 1, 4)

who_region <- read_csv("~/Desktop/GES/data/who-region.csv")
who_region <- who_region[,c("name", 'region')]

df <- merge(df, who_region, by.x='country', by.y='name', all.x=TRUE)
df$region <- ifelse(df$country=='Czech Republic', "Europe",
                    ifelse(df$country=="Russia", "Europe",
                           ifelse(df$country=="Taiwan", "Asia",
                                  ifelse(df$country=="United Kingdom", 'Europe',
                                         ifelse(df$country=="USA", "Americas", df$region))))) 

df <- merge(df, lengths, by.x='Contig', by.y='contig', all.x=TRUE)

df <- df %>% group_by(region, sampling.year) %>% mutate(region.value = n())
df <- df %>% group_by(genus, sampling.year) %>% mutate(genus.value = n())



library(ggplot2)
p.region <- ggplot(data=df, aes(x=sampling.year, fill=as.factor(region))) +
  geom_bar(stat="count", position='stack', size=0.5, color='black', width = 1) +
  labs(x = "", y='Count') +
  scale_fill_brewer("", palette="Set2") +
  scale_y_continuous(expand = c(0,0), limits=c(-0.05, 20)) +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=16),
        axis.ticks = element_line(size = 1, color='black'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color='black'))
  

p.genus <- ggplot(data=df, aes(x=sampling.year, fill=as.factor(genus))) +
  geom_bar(stat="count", position='stack', size=0.5, color='black', width = 1) +
  labs(x = "Sampling year", y='Count') +
  scale_fill_brewer("", palette="Set3") +
  scale_y_continuous(expand = c(0,0), limits=c(-0.05, 20)) +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=16),
        axis.ticks = element_line(size = 1, color='black'),
        axis.text.x = element_text(color='black', angle = 45, vjust = 0.8, hjust=0.8),
        axis.text.y = element_text(color='black'))

  

ggplot(data=df, aes(x=sampling.year, y=length)) +
  geom_boxplot()

summary(as.numeric(df$length))

summary(as.factor(df$region))
summary(as.factor(df$genus))

contig.count <- df %>% group_by(sampling.year) %>% summary(count = n())


library(cowplot)
grid::current.viewport()
cowplot::plot_grid(p.region, NULL, p.genus, rel_heights = c(1, -0.1, 1),
                   nrow=3, align='hv')

