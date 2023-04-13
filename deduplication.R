# first we read in the metadata table downloaded from ncbi
library(readr)
ncbiGES <- read_delim("~/Desktop/GES/data/NCBI-GES.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
# then we read in the list of contigs we confirmed to carry GES-5 using blastx
contigs <- read_csv("~/Desktop/GES/data/contigs-GES-5.txt", 
                                     col_names = FALSE)
# format the contig names
contigs$X1 <- gsub(".fasta", "", contigs$X1)

# next we filter the metadata for the GES-5 contigs
library(dplyr)
ncbiGES <- ncbiGES[ncbiGES$Contig %in% contigs$X1,]
# these contigs have multiple GES alleles (GES-5 + something else)
reps <- ncbiGES %>% group_by(Contig) %>% filter(n()>1)
# so we also have to drop the rows for the non-GES-5 alleles
# this gives us 1 row per contig, and each contig contains GES-5
ncbiGES <- ncbiGES[ncbiGES$`Element symbol`=="blaGES-5" ,]

# now we read in biosample/bioproject table
bioprojects <- read.csv("/Users/willmatlock/Desktop/GES/data/bioprojects.tsv", 
                          sep="\t", header = FALSE, fill=TRUE, na=NA)
# some wrangling is needed
bioprojects$V5 <- ""
bioprojects[c(34,36),"V5"] <- "PRJNA224116"
bioprojects <- bioprojects[-c(35,37),]
# we also need to manually curate bioprojects to remove meta-analyses labels
# PRJNA514245 - Pathogen Detection Assembly Project
# PRJNA224116 - Refseq Prokaryotic Genome Annotation Project
library(stringr)
library(tidyr)
bioprojects[bioprojects == ""] <- NA
bioprojects <- bioprojects %>%
  mutate_at(2:5, ~ replace(., str_detect(., 'PRJNA514245|PRJNA224116'), NA)) %>%
  # combine the bioproject labels into one string
  unite("bioproject", 2:5, sep = "-", na.rm=TRUE, remove = TRUE)
# we then combine our bioproject and ncbi tables
bioprojects <- merge(bioprojects, ncbiGES[c("BioSample", "Contig")], by.x='V1', by.y='BioSample', all.Y=TRUE)
bioprojects <- unique(bioprojects)
# rename columns
colnames(bioprojects) <- c("biosample", "bioproject", "contig")

# next we need to add contig lengths
lengths <- read_delim("~/Desktop/GES/data/lengths.tsv", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
colnames(lengths) <- c('contig', 'length')
# wrangle contig names from fasta headers
lengths$contig <- word(lengths$contig, 1)
# combine bioproject and lengths tables
bioprojects <- merge(bioprojects, lengths, by='contig', all.x=TRUE)
bioprojects <- unique(bioprojects)

# for sense checking after deduplication
no.biosamples <- length(unique(bioprojects$biosample))
no.bioprojects <- length(unique(bioprojects$bioproject))

# first we deduplicate by contig lengths
bioprojects <- bioprojects %>%
  # record the number of contigs in each bioproject
  group_by(bioproject) %>%
  mutate(bioproject.size = n()) %>%
  # then for each biosample
  group_by(biosample) %>%
  # we select longest contig
  # if multiple contigs are maximal, we choose a random representative
  mutate(longest = max(length),
         is.longest = ifelse(length==longest, TRUE, FALSE),
         biosample.size = n(),
         multiple.longest = ifelse(sum(is.longest)>1, TRUE, FALSE),
         rep = ifelse(multiple.longest==TRUE, sample(unique(contig), 1), contig[is.longest==TRUE])) %>%
  # filter output for our chosen representatives
  filter(contig==rep)

# sense check - we should still have the same number of bioprojects/biosamples as before we filtered
length(unique(bioprojects$biosample))==no.biosamples
length(unique(bioprojects$bioproject))==no.bioprojects

# second we deduplicate by sequence containment 
# all possible contig pairs
pairs <- merge(bioprojects$contig, bioprojects$contig) 
# rename columns
colnames(pairs) <- c("contig.x", "contig.y")
# add metadata to pair table
pairs <- merge(pairs, bioprojects, by.x='contig.x', by.y='contig', all.x=TRUE)
pairs <- merge(pairs, bioprojects, by.x='contig.y', by.y='contig', all.x=TRUE)
# filter down our pairs of representative contigs within the same bioproject
pairs <- pairs[pairs$bioproject.x == pairs$bioproject.y,] 

# now import mash screen table
mash <- read_delim("~/Desktop/GES/data/mash-screen.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          col_names = FALSE, trim_ws = TRUE)
mash <- mash[c(1,5)]
# rename columns
# NB the mash screen measures the containment of contig.x in contig.y
colnames(mash) <- c('x.in.y', 'contig.x')
# format contig names
mash$contig.x <- gsub('./contigs-GES-5/', '', mash$contig.x)
mash$contig.x <- gsub('.fasta', '', mash$contig.x)
# order of contigs used in mash screen
contigs.ordered <- mash[1:431,]$contig.x
mash$contig.y <- as.vector(sapply(contigs.ordered, function (x) rep(x,431)))

# combine mash and pairs data
containment <- merge(mash, pairs, all.y=TRUE)
# swap containment order for other half of data
colnames(mash) <- c('x.in.y', 'contig.y', 'contig.x')
containment <- merge(mash, pairs, all.y=TRUE)
containment <- drop_na(containment)

# first we collect contigs from bioprojects with one contig i.e., singletons
containment <- containment %>%
  # we can group by either bioproject.x or bioproject.y - they are the same now
  group_by(bioproject.x) %>%
  # number of contigs in each bioproject
  mutate(bioproject.n = n_distinct(c(contig.x, contig.y)),
         # does the bioproject contain only one contig
         is.singleton = ifelse(bioproject.n==1, TRUE, FALSE))
# our list of contigs to keep
# we can take either rep.x or rep.y - they are the same
to.keep <- containment[containment$is.singleton==TRUE,]$rep.x

# now we work on non-singleton bioprojects 
# if a bioproject has every contig contained in at least one other, we take the longest contig
# if multiple contigs are equally the longest, we choose a random representative 
# otherise, we can simply remove all the contigs which are contained in another
containment <- containment %>%
  # drop singletons and self-pairs
  filter(is.singleton==FALSE & rep.x!=rep.y) %>%
  # then within each bioproject
  group_by(bioproject.x) %>%
  # is contig.x fully contained in contig.y
  mutate(x.is.contained = ifelse(x.in.y==1, TRUE, FALSE),
         # how many contigs are contained in at least one other
         n.contained = n_distinct(contig.x[x.is.contained==TRUE]),
         # are all the contigs contained in at least one other
         all.contained = ifelse(bioproject.n == n.contained, TRUE, FALSE))
# so now we can list the contained contigs to drop  
to.remove <- unique(containment[containment$all.contained==FALSE & containment$x.in.y==1,]$rep.x)
# and remove them
containment <- containment %>%
  filter(!(rep.x %in% to.remove))
# leaving those to keep
to.keep <- unique(c(to.keep, containment[containment$all.contained==FALSE,]$rep.x))

# finally, let's work on the bioprojects where every contig is contained in at least one other
containment <- containment %>%
  filter(all.contained==TRUE) %>%
  group_by(bioproject.x) %>%
  # longest contig in each bioproject
  mutate(longest = max(length.y),
         is.longest = ifelse(length.y==longest, TRUE, FALSE)) %>%
  filter(is.longest==TRUE) %>%
  group_by(bioproject.x) %>%
  # sample a random longest contig
  mutate(rep = sample(unique(rep.y), 1))
# and add it to the rest
to.keep <- unique(c(to.keep, containment$rep))
  
# sense check - we should have the same number of bioprojects, but fewer biosamples
summary <- bioprojects[bioprojects$contig %in% to.keep,]
length(unique(summary$biosample))==no.biosamples
length(unique(summary$bioproject))==no.bioprojects
# we should also have one contig per biosample
length(unique(summary$contig))==length(unique(summary$biosample))

# write the list of deduplicated contigs with or without file extensions
write.csv(to.keep, '~/Desktop/contigs.txt', quote=FALSE, row.names=FALSE)
write.csv(paste0(to.keep, ".fasta"), '~/Desktop/contigs.txt', quote=FALSE, row.names=FALSE)




