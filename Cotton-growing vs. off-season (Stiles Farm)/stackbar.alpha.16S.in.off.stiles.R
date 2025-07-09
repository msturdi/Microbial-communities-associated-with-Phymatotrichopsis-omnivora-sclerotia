library(tidyverse)
library(readxl)
library(ggtext)
library(MoMAColors)
library(glue)
library(dplyr)

display.all.moma(24, colorblind_only=F)

metadata <- read_excel("metadata.16S.in.off.stiles.xlsx") %>%
  mutate(combo = glue("{season} {treatment}"))

otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.bacterial.stiles.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.bacterial.stiles.taxonomy") %>%
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
  mutate(season = factor(season, 
                         levels=c("in",
                                  "off")))

order_rel_abund <- otu_rel_abund %>%
  filter(level=="order") %>%
  group_by(season, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(season, taxon) %>%
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
  summarize(pool = max(rel_abund) < 0.019, #less than 2%
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
                    breaks=c("Caryophanales", "Gaiellales", "Rubrobacterales", "Hyphomicrobiales", "Solirubrobacterales", "Acidobacteria Gp6 incertae sedis", "Rhodospirillales", 
                             "Unclassified Actinomycetota", "Gemmatimonadales", "Micromonosporales", "Propionibacteriales", "Acidimicrobiales", "Micrococcales", "Unclassified Actinobacteria",
                             "Unclassified Thermoleophilia", "Myxococcales", "Sphingomonadales", "Burkholderiales", "Geodermatophilales", "Acidobacteria Gp16 incertae sedis",
                             "Streptomycetales", "Pseudonocardiales", "Unclassified Bacteria", "Other"),
                    labels=c("Caryophanales", "Gaiellales", "Rubrobacterales", "Hyphomicrobiales", "Solirubrobacterales", "Acidobacteria Gp6 incertae sedis", "Rhodospirillales", 
                             "Unclassified Actinomycetota", "Gemmatimonadales", "Micromonosporales", "Propionibacteriales", "Acidimicrobiales", "Micrococcales", "Unclassified Actinobacteria",
                             "Unclassified Thermoleophilia", "Myxococcales", "Sphingomonadales", "Burkholderiales", "Geodermatophilales", "Acidobacteria Gp16 incertae sedis",
                             "Streptomycetales", "Pseudonocardiales", "Unclassified Bacteria", "Other"),
                    values = c(moma.colors("Warhol", 23, direction=1, type="continuous"), "dimgrey")) +
  scale_x_discrete(limits=c("SI1A1B", "SI1A2B", "SI1B1B", "SI1B2B", "SI2A1B", "SI2A2B", "SI2B1B", "SI2B2B", "SI3A1B", "SI3A2B", "SI3B1B", "SI3B2B",
                            "SI1C1B", "SI1C2B", "SI2C1B", "SI2C2B", "SI3C1B", "SI3C2B",
                            "SIBA1B", "SIBA2B", "SIBB1B", "SIBB2B",
                            "SO1A1B", "SO1A2B", "SO1B1B", "SO1B2B", "SO2A1B", "SO2A2B", "SO2B1B", "SO2B2B", "SO3A1B", "SO3A2B", "SO3B1B", "SO3B2B",
                            "SO1C1B", "SO1C2B", "SO2C1B", "SO2C2B", "SO3C1B", "SO3C2B",
                            "SOBA1B", "SOBA2B", "SOBB1B", "SOBB2B"),
                   labels=c("In<br>Live<br>1", "In<br>Live<br>2", "In<br>Live<br>3", "In<br>Live<br>4", "In<br>Live<br>5", "In<br>Live<br>6", "In<br>Live<br>7", "In<br>Live<br>8", 
                            "In<br>Live<br>9", "In<br>Live<br>10", "In<br>Live<br>11", "In<br>Live<br>12",
                            "In<br>Heat<br>killed<br>1", "In<br>Heat<br>killed<br>2", "In<br>Heat<br>killed<br>3", "In<br>Heat<br>killed<br>4", "In<br>Heat<br>killed<br>5", 
                            "In<br>Heat<br>killed<br>6",
                            "In<br>Bulk<br>1", "In<br>Bulk<br>2", "In<br>Bulk<br>3", "In<br>Bulk<br>4",
                            "Off<br>Live<br>1", "Off<br>Live<br>2", "Off<br>Live<br>3", "Off<br>Live<br>4", "Off<br>Live<br>5", "Off<br>Live<br>6", "Off<br>Live<br>7", "Off<br>Live<br>8", 
                            "Off<br>Live<br>9", "Off<br>Live<br>10", "Off<br>Live<br>11", "Off<br>Live<br>12",
                            "Off<br>Heat<br>killed<br>1", "Off<br>Heat<br>killed<br>2", "Off<br>Heat<br>killed<br>3", "Off<br>Heat<br>killed<br>4", "Off<br>Heat<br>killed<br>5", 
                            "Off<br>Heat<br>killed<br>6",
                            "Off<br>Bulk<br>1", "Off<br>Bulk<br>2", "Off<br>Bulk<br>3", "Off<br>Bulk<br>4")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL,
       y="Relative Abundance (%)",
       title="Order-level Bacterial Alpha Diversity: Off-season vs. 2024 Cotton-growing Season<br>**Stiles Farm**") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.16S.in.off.stiles.png", width=22, height=9)

#order-level averages
inner_join(order_rel_abund_mean, order_pool_mean, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(combo, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), 
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ggplot(aes(x=combo, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name="Orders",
                    breaks=c("Gaiellales", "Caryophanales", "Rubrobacterales", "Hyphomicrobiales", "Solirubrobacterales", "Acidobacteria Gp6 incertae sedis", "Rhodospirillales", 
                             "Unclassified Actinomycetota", "Gemmatimonadales", "Micromonosporales", "Propionibacteriales", "Acidimicrobiales", "Micrococcales", "Unclassified Actinobacteria",
                             "Unclassified Thermoleophilia", "Sphingomonadales", "Myxococcales", "Burkholderiales", "Geodermatophilales", "Streptomycetales", "Acidobacteria Gp16 incertae sedis", 
                             "Pseudonocardiales", "Acidobacteria Gp3 incertae sedis", "Unclassified Bacteria", "Other"),
                    labels=c("**Gaiellales (Off)***", "**Caryophanales (In)***", "Rubrobacterales", "**Hyphomicrobiales (Off)***", "Solirubrobacterales", 
                             "**Acidobacteria Gp6 incertae sedis (Off)***", "Rhodospirillales", "**Unclassified Actinomycetota (In)***", "**Gemmatimonadales (In)***", "Micromonosporales", 
                             "Propionibacteriales", "Acidimicrobiales", "Micrococcales", "Unclassified Actinobacteria", "Unclassified Thermoleophilia", "Sphingomonadales", "Myxococcales", 
                             "**Burkholderiales (Off)***", "Geodermatophilales", "**Streptomycetales (In)***", "Acidobacteria Gp16 incertae sedis", "**Pseudonocardiales (In)***", 
                             "**Acidobacteria Gp3 incertae sedis (In)***", "**Unclassified Bacteria (In)***", "Other"),
                    values = c(moma.colors("Warhol", 24, direction=1, type="continuous"), "dimgrey")) +
  scale_x_discrete(limits=c("in live", "in dead", "in none", "off live", "off dead", "off none"),
                   labels=c("In-season<br>Live (n=12)", "In-season<br>Heat-killed<br>(n=6)", "In-season<br>Bulk (n=4)",
                            "Off-season<br>Live (n=12)", "Off-season<br>Heat-killed<br>(n=6)", "Off-season<br>Bulk (n=4)")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL,
       y="Average Relative Abundance (%)",
       title="Order-level Bacterial Alpha Diversity<br>In- and Off-season Stiles Farm") +
  labs(caption = "Other = relative abundance < 1%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        plot.caption = element_text(hjust=0.95, size=8),
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.16S.in.off.stiles.ave.png", width=13, height=8)
