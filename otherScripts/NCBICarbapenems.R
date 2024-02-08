# Filtering MicroBIGG-E table 10/2/2023 using BigQuery using SQL command:
#
# SELECT
# taxgroup_name,
# element_symbol,
# LEFT(collection_date, 4),
# COUNT(*),
# FROM
# `ncbi-pathogen-detect.pdbrowser.microbigge`
# WHERE
# subclass LIKE '%CARBAPENEM%'
# AND subtype = 'AMR'
# AND LEFT(collection_date, 4) IS NOT NULL
# GROUP BY
# taxgroup_name,
# element_symbol,
# LEFT(collection_date, 4)
#

library(readr)
data <- read_csv("~/Desktop/ncbi-carbapenems.csv")
colnames(data) <- c('taxgroup', 'gene', 'year', 'count')

# merge alleles and clean genera
library(dplyr)
data <- data %>%
  mutate(gene = gsub("-.*", "", gene)) %>%
  mutate(taxgroup = gsub('E.coli and Shigella', 'Escherichia/Shigella', taxgroup)) %>%
  mutate(taxgroup = gsub('Kluyvera_intermedia', 'Kluyvera', taxgroup)) %>%
  mutate(taxgroup = gsub(" .*", "", taxgroup)) %>%
  group_by(taxgroup, gene, year) %>%
  mutate(count = sum(count)) %>%
  distinct()

data <- data[data$year>=1990,]

# plot timeline
library(ggplot2)
plot <- ggplot(data) + 
  geom_jitter(aes(x = year, y = gene, colour = as.factor(taxgroup), size=count), alpha=0.8, width=0.1)
plot


