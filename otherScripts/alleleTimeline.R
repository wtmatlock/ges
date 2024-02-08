library(readr)

alleles <- read_delim("~/Desktop/blaGES/data//allele_counts_by_year.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

alleles <- alleles[alleles$Ambler_class=="bla-A",] # filter to class A bla
#alleles <- alleles[alleles$`2022`>=20,] # filter to >= 10 in 2022
alleles <- alleles[alleles$Subclass!='BETA-LACTAM',] # filter out non-ESBL/carbapenemase

library(reshape2)

alleles <- melt(alleles)
alleles$blaGES <- ifelse(alleles$Family=="blaGES", TRUE, FALSE)
alleles$Label <- ifelse(alleles$variable==2022, as.character(alleles$Family), NA_character_)

alleles$Subclass <- ifelse(alleles$Subclass=="CARBAPENEM", "Carbapenemase", "ESBL")

colnames(alleles) <- c("Gene.family", "Hydrolytic.profile", "Ambler.class", "Year", "Count", "blaGES", "Label")

library(dplyr)

alleles <- alleles %>%
  group_by(Gene.family, Hydrolytic.profile, Ambler.class, blaGES) %>%
  mutate(Cumulative.count = cumsum(Count))

library(ggplot2)
library(ggrepel)

p <- ggplot(data=alleles, aes(x=Year, y=Count, group=Gene.family, colour=blaGES)) +
  scale_color_manual(values=c('#377eb8', '#e41a1c')) +
  geom_line(alpha=0.5, linewidth=1) +
  theme_linedraw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(Hydrolytic.profile), scales="free") +
  geom_label_repel(aes(label = Label),
                   nudge_x = 1,
                   na.rm = TRUE)
p

