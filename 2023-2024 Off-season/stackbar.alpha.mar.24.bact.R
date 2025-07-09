library(tidyverse)
library(readxl)
library(ggtext)
library(MoMAColors)
library(glue)
library(dplyr)

display.all.moma(24, colorblind_only=F)

metadata <- read_excel("mar.24.bact.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}"))

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.mar.24.bact.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.mar.24.bact.taxonomy") %>%
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
                                 "_", " "))

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
  summarize(pool = max(mean_rel_abund) < 1,
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
                    breaks=c("Gaiellales", "Unclassified Bacteria", "Caryophanales", "Acidobacteria Gp6 incertae sedis", "Solirubrobacterales", "Rhodospirillales",
                             "Rubrobacterales", "Micromonosporales", "Propionibacteriales", "Micrococcales", "Unclassified Actinomycetota", "Acidimicrobiales", "Myxococcales",
                             "Unclassified Thermoleophilia", "Gemmatimonadales", "Sphingomonadales", "Burkholderiales", "Acidobacteria Gp16 incertae sedis",
                             "Unclassified Actinobacteria", "Geodermatophilales", "Streptomycetales", "Pseudonocardiales", "Chitinophagales", "Acidobacteria Gp4 incertae sedis",
                             "Nevskiales", "Hyphomicrobiales", "Other"),
                    labels=c("Gaiellales", "Unclassified Bacteria", "Caryophanales", "Acidobacteria Gp6 incertae sedis", "Solirubrobacterales", "Rhodospirillales",
                             "Rubrobacterales", "Micromonosporales", "Propionibacteriales", "Micrococcales", "Unclassified Actinomycetota", "Acidimicrobiales", "Myxococcales",
                             "Unclassified Thermoleophilia", "Gemmatimonadales", "Sphingomonadales", "Burkholderiales", "Acidobacteria Gp16 incertae sedis",
                             "Unclassified Actinobacteria", "Geodermatophilales", "Streptomycetales", "Pseudonocardiales", "Chitinophagales", "Acidobacteria Gp4 incertae sedis",
                             "Nevskiales", "Hyphomicrobiales", "Other"),
                    values = c(moma.colors("Warhol", 26, direction=1, type="continuous"), "dimgrey")) +
  scale_x_discrete(limits=c("BFO1A1B", "BFO1A2B", "BFO1B1B", "BFO1B2B", "BFO2A1B", "BFO2A2B", "BFO2B1B", "BFO2B2B", "BFO3A1B", "BFO3A2B", "BFO3B1B", "BFO3B2B",
                            "BFO1C1B", "BFO1C2B", "BFO2C1B", "BFO2C2B", "BFO3C1B", "BFO3C2B",
                            "BFOBA1B", "BFOBA2B", "BFOBB1B", "BFOBB2B",
                            "SO1A1B", "SO1A2B", "SO1B1B", "SO1B2B", "SO2A1B", "SO2A2B", "SO2B1B", "SO2B2B", "SO3A1B", "SO3A2B", "SO3B1B", "SO3B2B",
                            "SO1C1B", "SO1C2B", "SO2C1B", "SO2C2B", "SO3C1B", "SO3C2B",
                            "SOBA1B", "SOBA2B", "SOBB1B", "SOBB2B"),
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
       title="Order-level Bacterial Alpha Diversity: **2023-2024 Off-season**") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.16S.png", width=18, height=8)
