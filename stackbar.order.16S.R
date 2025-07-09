library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(MoMAColors)
library(glue)
library(dplyr)

display.all.moma(16, colorblind_only=F)

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

order_rel_abund <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(trt, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(trt, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " ")) #gets ride of _ in taxon names

order_pool <- order_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(rel_abund) < 0.025, 
            mean = mean(rel_abund),
            .groups="drop")

inner_join(order_rel_abund, order_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), 
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>% #all of this mutate step is creating an "anchor" at the top and the bottom of the figure using the 2 most abundant phyla
  ggplot(aes(x=sample, y=rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name="Orders",
                    breaks=c("Unclassified Actinobacteria", "Rubrobacterales", "Unclassified Bacteria", "Rhizobiales", "Bacillales",
                             "Acidobacteria Gp6 order incertae sedis", "Solirubrobacterales", "Rhodospirillales", "Acidimicrobiales",
                             "Alphaproteobacteria order incertae sedis", "Sphingomonadales", "Burkholderiales", "Sphingobacteriales", 
                             "Pseudomonadales", "Actinomycetales", "Other"),
                    labels=c("Unclassified Actinobacteria", "Rubrobacterales", "Unclassified Bacteria", "Rhizobiales", "Bacillales",
                             "Acidobacteria Gp6 order incertae sedis", "Solirubrobacterales", "Rhodospirillales", "Acidimicrobiales",
                             "Alphaproteobacteria order incertae sedis", "Sphingomonadales", "**Burkholderiales***", "Sphingobacteriales", 
                             "Pseudomonadales", "Actinomycetales", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Warhol", 15, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
  scale_x_discrete(limits=c("AA1B", "AB1B", "AA2B", "AB2B", "AA3B", "AB3B", "BA1B", "BB1B", "BA2B", "BB2B", "BA3B", "BB3B", 
                            "CA1B", "CB1B", "CA2B", "CB2B", "CA3B", "CB3B", 
                            "BAAB", "BABB", "BBAB", "BBBB"),
                   labels=c("Live<br>1", "Live<br>2", "Live<br>3", "Live<br>4", "Live<br>5","Live<br>6", "Live<br>7", "Live<br>8", "Live<br>9", "Live<br>10",
                            "Live<br>11", "Live<br>12",
                            "Heat<br>Killed<br>1", "Heat<br>Killed<br>2", "Heat<br>Killed<br>3", "Heat<br>Killed<br>4", "Heat<br>Killed<br>5", 
                            "Heat<br>Killed<br>6",
                            "Bulk<br>Soil<br>1", "Bulk<br>Soil<br>2", "Bulk<br>Soil<br>3", "Bulk<br>Soil<br>4")) +
  scale_y_continuous(expand=c(0,0)) + #this makes the bars touch the y axis
  labs(x=NULL,
       y="Relative Abundance (%)",
       title="Order-level 16S rRNA Gene Alpha Diversity") +
  labs(caption = "Other = relative abundance < 2.5%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        plot.caption = element_text(hjust=1.3, size=8),
        plot.title = element_text(hjust=0.5))

ggsave("order.stackedbar.16S.png", width=12, height=5)

#class-level
class_rel_abund <- otu_rel_abund %>%
  filter(level=="class") %>%
  group_by(trt, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(trt, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " ")) #gets ride of _ in taxon names

class_pool <- class_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(rel_abund) < 0.02, 
            mean = mean(rel_abund),
            .groups="drop")

inner_join(class_rel_abund, class_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), 
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>% #all of this mutate step is creating an "anchor" at the top and the bottom of the figure using the 2 most abundant phyla
  ggplot(aes(x=sample, y=rel_abund, fill=taxon)) +
  geom_col() 
