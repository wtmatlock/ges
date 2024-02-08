library(ggtree)
library(phangorn)

tree <- read.tree("/Users/willmatlock/Desktop/GES/ges-allele-phylo/NCBI-GES.faa.phb")
tree <- midpoint(tree)

hydrolytic <- as.data.frame(tree$tip.label)
carb <- c("GES-2", "GES-50", "GES-24", "GES-48",
          "GES-47", "GES-54", "GES-4", "GES-51",
          "GES-49", "GES-39", "GES-6", "GES-53",
          "GES-15", "GES-14", "GES-20", "GES-16",
          "GES-5", "GES-21", "GES-18")
esbl <- c("GES-1", "GES-8", "GES-23", "GES-12",
          "GES-9", "GES-11", "GES-22", "GES-10",
          "GES-3", "GES-13", "GES-7", "GES-17",
          "GES-52", "GES-38", "GES-19")
hydrolytic$func <- ifelse(hydrolytic$`tree$tip.label` %in% carb, "Carbapenemase",
                          ifelse(hydrolytic$`tree$tip.label` %in% esbl, "ESBL", "Unknown"))

### alignment tree

library(readr)
library(dplyr)
library(stringr)

df <- read_delim("~/Desktop/GES/data/NCBI-GES.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

summary(as.factor(df$`Element symbol`))

df$genera <- word(df$`#Scientific name`, 1)
df$country <- word(df$Location, sep = "\\:")
df$year <- str_sub(df$`Collection date`, 1, 4)   

df <- df %>%
  group_by(`Element symbol`, year, country, genera) %>%
  summarise(count = n())

df$`Element symbol` <- gsub("bla", "", df$`Element symbol`)

df <- merge(df, hydrolytic, by.x='Element symbol', by.y='tree$tip.label', all.y=TRUE)

### tree

p.tree <- ggtree(tree) + geom_tiplab(align=TRUE) + theme_tree2() + ggplot2::xlim(0, 0.04)

df$element_order <- factor(df$`Element symbol`, levels = rev(get_taxa_name(p1))

### genera

library(ggplot2)
library(ggstance)
library(ggnewscale)
library(ggh4x)

p.hydrolytic.profile <- ggplot(df, aes(x=element_order)) + 
  geom_bar(aes(y=..count.., fill=func), position='fill',size=0.25, color='black', width = 1) + 
  coord_flip() +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  labs(x = "", y='') +
  scale_fill_manual("", values = c("#de2d26", "#fcbba1", "white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8)) + 
  theme(legend.position = "none")

df2 <- df %>% group_by(element_order, genera) %>% summarise(sum = sum(count))

p.genus <- ggplot(df2, aes(x=element_order, y=sum, fill=genera)) + 
  geom_bar(stat="identity" ,size=0.25, color='black', width = 1) + 
  coord_flip() +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400)) +
  labs(x = "", y='') +
  scale_fill_brewer("", palette = 'Set3') +
  theme(panel.grid.major = element_line(colour = "grey"),
        panel.grid.minor = element_line(colour = "grey"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=0.8))

library(cowplot)
grid::current.viewport()
cowplot::plot_grid(p.tree,
                   p.hydrolytic.profile,
                   p.genus,
                   nrow=1, align='h',
                   rel_widths = c(3, 1, 5),
                   scale=c(1, 1.01, 1.01),
                   labels = c("a", "b", "c"))

