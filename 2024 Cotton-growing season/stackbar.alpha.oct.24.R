library(tidyverse)
library(readxl)
library(ggtext)
library(MoMAColors)
library(glue)
library(dplyr)

display.all.moma(24, colorblind_only=F)

metadata <- read_excel("oct.24.16S.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}"))

shared <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.bacterial.shared",
                   col_types = cols(Group = col_character(),
                                    .default = col_double()))

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.bacterial.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.oct.24.bacterial.taxonomy") %>%
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
  mutate(location = factor(location, 
                           levels=c("BF",
                                    "Stiles")))

order_rel_abund <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(location, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(location, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " ")) #gets ride of _ in taxon names

order_rel_abund_mean <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(combo, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " "))

order_pool <- order_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(rel_abund) < 0.019, #less than 1.9%
            mean = mean(rel_abund),
            .groups="drop")

order_pool_mean <- order_rel_abund_mean %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, #less than 1% - this is not in decimal percentage because with the means I multiplied by 100
            mean = mean(mean_rel_abund),
            .groups="drop")

#order-level
inner_join(order_rel_abund, order_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), 
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ggplot(aes(x=sample, y=rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name="Orders",
                    breaks=c("Caryophanales", "Hyphomicrobiales", "Rubrobacterales", "Solirubrobacterales", "Gaiellales", "Rhodospirillales", "Acidobacteria Gp6 incertae sedis",
                             "Unclassified Actinomycetota", "Micromonosporales", "Propionibacteriales", "Gemmatimonadales", "Acidimicrobiales", "Unclassified Thermoleophilia",
                             "Micrococcales", "Myxococcales", "Pseudonocardiales", "Acidobacteria Gp16 incertae sedis", "Streptomycetales", "Unclassified Actinobacteria",
                             "Sphingomonadales", "Acidobacteria Gp4 incertae sedis", "Streptosporangiales", "Nevskiales", "Coleofasciculales", "Unclassified Bacteria",
                             "Other"),
                    labels=c("Caryophanales", "Hyphomicrobiales", "Rubrobacterales", "Solirubrobacterales", "Gaiellales", "Rhodospirillales", 
                             "Acidobacteria Gp6 incertae sedis", "Unclassified Actinomycetota", "Micromonosporales", "Propionibacteriales", 
                             "Gemmatimonadales", "Acidimicrobiales", "Unclassified Thermoleophilia", "Micrococcales", "Myxococcales", 
                             "Pseudonocardiales", "Acidobacteria Gp16 incertae sedis", "Streptomycetales", "Unclassified Actinobacteria",
                             "Sphingomonadales", "Acidobacteria Gp4 incertae sedis", "Streptosporangiales", "Nevskiales", 
                             "Coleofasciculales", "Unclassified Bacteria", "Other"),
                    values = c(moma.colors("Warhol", 25, direction=1, type="continuous"), "dimgrey")) +
  scale_x_discrete(limits=c("BI1A1B", "BI1A2B", "BI1B1B", "BI1B2B", "BI2A1B", "BI2A2B", "BI2B1B", "BI2B2B", "BI3A1B", "BI3A2B", "BI3B1B", "BI3B2B",
                            "BI1C1B", "BI1C2B", "BI2C1B", "BI2C2B", "BI3C1B", "BI3C2B",
                            "BIBA1B", "BIBA2B", "BIBB1B", "BIBB2B",
                            "SI1A1B", "SI1A2B", "SI1B1B", "SI1B2B", "SI2A1B", "SI2A2B", "SI2B1B", "SI2B2B", "SI3A1B", "SI3A2B", "SI3B1B", "SI3B2B",
                            "SI1C1B", "SI1C2B", "SI2C1B", "SI2C2B", "SI3C1B", "SI3C2B",
                            "SIBA1B", "SIBA2B", "SIBB1B", "SIBB2B"),
                   labels=c("BF<br>Live<br>1", "BF<br>Live<br>2", "BF<br>Live<br>3", "BF<br>Live<br>4", "BF<br>Live<br>5","BF<br>Live<br>6", "BF<br>Live<br>7", "BF<br>Live<br>8", 
                            "BF<br>Live<br>9", "BF<br>Live<br>10", "BF<br>Live<br>11", "BF<br>Live<br>12",
                            "BF<br>Heat<br>Killed<br>1", "BF<br>Heat<br>Killed<br>2", "BF<br>Heat<br>Killed<br>3", "BF<br>Heat<br>Killed<br>4", "BF<br>Heat<br>Killed<br>5", 
                            "BF<br>Heat<br>Killed<br>6",
                            "BF<br>Bulk<br>Soil<br>1", "BF<br>Bulk<br>Soil<br>2", "BF<br>Bulk<br>Soil<br>3", "BF<br>Bulk<br>Soil<br>4",
                            "SF<br>Live<br>1", "SF<br>Live<br>2", "SF<br>Live<br>3", "SF<br>Live<br>4", "SF<br>Live<br>5","SF<br>Live<br>6", "SF<br>Live<br>7", "SF<br>Live<br>8", 
                            "SF<br>Live<br>9", "SF<br>Live<br>10", "SF<br>Live<br>11", "SF<br>Live<br>12",
                            "SF<br>Heat<br>Killed<br>1", "SF<br>Heat<br>Killed<br>2", "SF<br>Heat<br>Killed<br>3", "SF<br>Heat<br>Killed<br>4", "SF<br>Heat<br>Killed<br>5", 
                            "SF<br>Heat<br>Killed<br>6",
                            "SF<br>Bulk<br>Soil<br>1", "SF<br>Bulk<br>Soil<br>2", "SF<br>Bulk<br>Soil<br>3", "SF<br>Bulk<br>Soil<br>4")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL,
       y="Relative Abundance (%)",
       title="Order-level Bacterial Alpha Diversity: **2024 Cotton-growing season**") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.16S.png", width=18, height=8)
