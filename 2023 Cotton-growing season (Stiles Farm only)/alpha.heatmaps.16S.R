library(tidyverse)
library(readxl)
library(ggtext)
library(glue)
library(dplyr)

metadata <- read_excel("16S.2023.metadata.xlsx")

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern='["]', replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), 
           sep=";")

otu_rel_abund <- inner_join(metadata, otu_counts, by="sample") %>%
  inner_join(., taxonomy, by="otu") %>%
  group_by(sample) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to="level",
               values_to="taxon") %>%
  mutate(trt = factor(trt, 
                      levels=c("live",
                               "dead",
                               "none")))

#order-level heat map
order_rel_abund <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(trt, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(trt, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"),
         taxon = str_replace_all(taxon, 
                                 "_", " "))

order_rel_abund %>%
  group_by(taxon) %>%
  summarize(rel_abund=sum(rel_abund)) %>%
  arrange(rel_abund) %>%
  view()

top_order <- order_rel_abund %>%
  group_by(taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  arrange(rel_abund) %>%
  top_n(14, rel_abund) %>%
  pull(taxon)
  view(top_order)

order_rel_abund %>%
  filter(taxon %in% top_order) %>%
  mutate(taxon=factor(taxon, levels=top_order)) %>%
  mutate(rel_abund= 100 * (rel_abund + 1/20000))

order_rel_abund %>%
  filter(taxon %in% top_order) %>%
  mutate(taxon=factor(taxon, levels=top_order)) %>%
  mutate(rel_abund= 100 * (rel_abund + 1/20000)) %>%
  ggplot(aes(x=sample, y=taxon, fill=rel_abund)) +
  geom_tile(color="black") +
  scale_fill_gradient(low="white", high="red",
                      name=("Relative \nAbundance (%)"),
                      limits=c(0,NA)) +
  coord_fixed() +
  scale_x_discrete(limits=c("AA1B", "AB1B", "AA2B", "AB2B", "AA3B", "AB3B", "BA1B", "BB1B", "BA2B", "BB2B", "BA3B", "BB3B", 
                            "CA1B", "CB1B", "CA2B", "CB2B", "CA3B", "CB3B", 
                            "BAAB", "BABB", "BBAB", "BBBB"),
                   labels=c("Live<br>1", "Live<br>2", "Live<br>3", "Live<br>4", "Live<br>5","Live<br>6", "Live<br>7", "Live<br>8", "Live<br>9", "Live<br>10",
                            "Live<br>11", "Live<br>12",
                            "Heat<br>Killed<br>1", "Heat<br>Killed<br>2", "Heat<br>Killed<br>3", "Heat<br>Killed<br>4", "Heat<br>Killed<br>5", 
                            "Heat<br>Killed<br>6",
                            "Bulk<br>Soil<br>1", "Bulk<br>Soil<br>2", "Bulk<br>Soil<br>3", "Bulk<br>Soil<br>4")) +
  labs(title="Relative Abundance of the Top 14 Orders",
       x=NULL,
       y=NULL) +
  theme_classic() +
  theme(axis.line = element_blank(), #gets rid of the axis lines
        axis.ticks = element_blank(), #gets rid of the axis tick marks
        axis.text.y = element_markdown(size=8),
        axis.text.x = element_markdown(size=8),
        #legend.title.align = 0.5,
        legend.title = element_text(size=8),
        legend.text = element_text(size=7.5),
        legend.key.height = unit(11, "pt"),
        plot.title = element_text(hjust=0.5))

ggsave("heatmap.16S.order.png", width=11, height=6)
