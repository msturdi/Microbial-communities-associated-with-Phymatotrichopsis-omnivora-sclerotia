library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(MoMAColors)
library(glue)
library(dplyr)

display.all.moma(24, colorblind_only=F)

metadata <- read_excel("oct.24.16S.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}"))

shared <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.bacterial.shared",
                   col_types = cols(Group = col_character(),
                                    .default = col_double()))

meta_cols <- c("label", "Group", "numOtus") #define columns to keep from shared file
otu_sums <- colSums(shared[, !(names(shared) %in% meta_cols)]) #sum OTU abundances across all samples
otus_to_keep <- names(otu_sums[otu_sums > 1]) #keep OTUs with total counts above 1

shared_filtered <- shared %>%
  select(all_of(meta_cols), all_of(otus_to_keep))

otu_counts <- shared_filtered %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count")

#otu_counts <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.bacterial.shared") %>%
  #select(Group, starts_with("Otu")) %>%
  #rename(sample = Group) %>%
  #pivot_longer(-sample, names_to="otu", values_to = "count")

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
         taxon = fct_shift(taxon, n=1)) %>% #all of this mutate step is creating an "anchor" at the top and the bottom of the figure using the 2 most abundant phyla
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
                    #labels=c("**Caryophanales (SF)***", "Hyphomicrobiales", "**Rubrobacterales (SF)***", "Solirubrobacterales", "**Gaiellales (SF)***", "**Rhodospirillales (BF)***", 
                             #"**Acidobacteria Gp6 incertae sedis (BF)***", "**Unclassified Actinomycetota (SF)***", "Micromonosporales", "**Propionibacteriales (BF)***", 
                             #"**Gemmatimonadales (SF)***", "**Acidimicrobiales (BF)***", "**Unclassified Thermoleophilia (BF)***", "**Micrococcales (BF)***", "**Myxococcales (BF)***", 
                             #"**Pseudonocardiales (BF)***", "**Acidobacteria Gp16 incertae sedis (BF)***", "**Streptomycetales (BF)***", "**Unclassified Actinobacteria (SF)***",
                             #"**Sphingomonadales (SF)***", "**Acidobacteria Gp4 incertae sedis (BF)***", "**Streptosporangiales (BF)***", "**Nevskiales (BF)***", 
                             #"**Coleofasciculales (BF)***", "**Unclassified Bacteria (SF)***", "Other"),
                    values = c(moma.colors("Warhol", 25, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
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
  #labs(caption = "Other = relative abundance < 2%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        #plot.caption = element_text(hjust=0.8, size=8),
        plot.title = element_markdown(hjust=0.5))

ggsave("order.stackedbar.16S.png", width=18, height=8)

#order-level averages
inner_join(order_rel_abund_mean, order_pool_mean, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(combo, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), 
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>% #all of this mutate step is creating an "anchor" at the top and the bottom of the figure using the 2 most abundant phyla
  ggplot(aes(x=combo, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name="Orders",
                    breaks=c("Caryophanales", "Rubrobacterales", "Hyphomicrobiales", "Solirubrobacterales", "Gaiellales", "Rhodospirillales", "Acidobacteria Gp6 incertae sedis",
                             "Unclassified Actinomycetota", "Micromonosporales", "Propionibacteriales", "Acidimicrobiales", "Gemmatimonadales", "Unclassified Thermoleophilia",
                             "Micrococcales", "Myxococcales", "Acidobacteria Gp16 incertae sedis", "Pseudonocardiales", "Streptomycetales", "Unclassified Actinobacteria",
                             "Geodermatophilales", "Sphingomonadales", "Burkholderiales", "Acidobacteria Gp4 incertae sedis", "Acidobacteria Gp3 incertae sedis", "Streptosporangiales", 
                             "Unclassified Bacteria", "Other"),
                    labels=c("**Caryophanales (SF)***", "**Rubrobacterales (SF)***", "Hyphomicrobiales", "Solirubrobacterales", "**Gaiellales (SF)***", "**Rhodospirillales (BF)***", 
                             "**Acidobacteria Gp6 incertae sedis (BF)***", "**Unclassified Actinomycetota (SF)***", "Micromonosporales", "**Propionibacteriales (BF)***", 
                             "**Acidimicrobiales (BF)***", "**Gemmatimonadales (SF)***", "**Unclassified Thermoleophilia (BF)***", "**Micrococcales (BF)***", "**Myxococcales (BF)***", 
                             "**Acidobacteria Gp16 incertae sedis (BF)***", "**Pseudonocardiales (BF)***", "**Streptomycetales (BF)***", "**Unclassified Actinobacteria (SF)***",
                             "Geodermatophilales", "**Sphingomonadales (SF)***", "Burkholderiales", "**Acidobacteria Gp4 incertae sedis (BF)***", "**Acidobacteria Gp3 incertae sedis (SF)***", 
                             "**Streptosporangiales (BF)***", "**Unclassified Bacteria (SF)***", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Warhol", 26, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
  scale_x_discrete(limits=c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                   labels=c("Bottom Farm<br>Live (n=12)", "Bottom Farm<br>Heat-killed<br>(n=6)", "Bottom Farm<br>Bulk (n=4)",
                            "Stiles Farm<br>Live (n=12)", "Stiles Farm<br>Heat-killed<br>(n=6)", "Stiles Farm<br>Bulk (n=4)")) +
  scale_y_continuous(expand=c(0,0)) + #this makes the bars touch the y axis
  labs(x=NULL,
       y="Average Relative Abundance (%)",
       title="Order-level Bacterial Alpha Diversity") +
  labs(caption = "Other = relative abundance < 1%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.position = "bottom",
        plot.caption = element_text(hjust=0.95, size=8),
        plot.title = element_text(hjust=0.5))

ggsave("order.stackedbar.16S.ave.png", width=13, height=8)

#class-level
class_rel_abund <- otu_rel_abund %>%
  filter(level=="class") %>%
  group_by(location, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(location, taxon) %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " ")) #gets ride of _ in taxon names

class_rel_abund_mean <- otu_rel_abund %>%
  filter(level=="class") %>%
  group_by(combo, sample, taxon) %>%
  summarize(rel_abund = sum(rel_abund), 
            .groups="drop") %>%
  group_by(combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified \\1"),
         taxon = str_replace_all(taxon, 
                                 "_", " "))

class_pool <- class_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(rel_abund) < 0.02, 
            mean = mean(rel_abund),
            .groups="drop")

class_pool_mean <- class_rel_abund_mean %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, #less than 1% - this is not in decimal percentage because with the means I multiplied by 100
            mean = mean(mean_rel_abund),
            .groups="drop")

#showing all samples
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
  geom_col() +
  scale_fill_manual(name="Classes",
                    breaks=c("Thermoleophilia", "Actinobacteria", "Unclassified Bacteria", "Bacilli", "Acidobacteria Gp6", "Rubrobacteria",
                             "Deltaproteobacteria", "Unclassified Actinomycetota", "Betaproteobacteria", "Acidimicrobiia", "Gammaproteobacteria",
                             "Gemmatimonadia", "Acidobacteria Gp16", "Chitinophagia", "Acidobacteria Gp4", "Alphaproteobacteria", "Other"),
                    labels=c("Thermoleophilia", "Actinobacteria", "Unclassified Bacteria", "Bacilli", "Acidobacteria Gp6", "Rubrobacteria",
                             "Deltaproteobacteria", "Unclassified Actinomycetota", "Betaproteobacteria", "Acidimicrobiia", "Gammaproteobacteria",
                             "Gemmatimonadia", "Acidobacteria Gp16", "Chitinophagia", "Acidobacteria Gp4", "Alphaproteobacteria", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Warhol", 16, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
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
       title="Class-level 16S rRNA Gene Alpha Diversity") +
  labs(caption = "Other = relative abundance < 2%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        plot.caption = element_text(hjust=1.14, size=8),
        plot.title = element_text(hjust=0.5))

ggsave("class.stackedbar.16S.png", width=18, height=8)

#showing means by combo
inner_join(class_rel_abund_mean, class_pool_mean, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(combo, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund), 
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>% #all of this mutate step is creating an "anchor" at the top and the bottom of the figure using the 2 most abundant phyla
  ggplot(aes(x=combo, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name="Classes",
                    breaks=c("Alphaproteobacteria", "Actinobacteria", "Unclassified Bacteria", "Bacilli", "Acidobacteria Gp6", "Rubrobacteria",
                             "Unclassified Actinomycetota", "Deltaproteobacteria", "Betaproteobacteria", "Acidimicrobiia", "Gemmatimonadia",
                             "Gammaproteobacteria", "Acidobacteria Gp16", "Chitinophagia", "Acidobacteria Gp4", "Thermoleophilia", "Other"),
                    labels=c("Alphaproteobacteria", "Actinobacteria", "Unclassified Bacteria", "Bacilli", "Acidobacteria Gp6", "Rubrobacteria",
                             "Unclassified Actinomycetota", "Deltaproteobacteria", "Betaproteobacteria", "Acidimicrobiia", "Gemmatimonadia",
                             "Gammaproteobacteria", "Acidobacteria Gp16", "Chitinophagia", "Acidobacteria Gp4", "Thermoleophilia", "Other"),
                    #values = c(brewer.pal(12, "Paired"), "navyblue", "dimgrey")) +
                    values = c(moma.colors("Warhol", 16, direction=1, type="continuous"), "dimgrey")) + #type could be continuous
  scale_x_discrete(limits=c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                   labels=c("Bottom Farm<br>Live (n=12)", "Bottom Farm<br>Heat-killed<br>(n=6)", "Bottom Farm<br>Bulk (n=4)",
                            "Stiles Farm<br>Live (n=12)", "Stiles Farm<br>Heat-killed<br>(n=6)", "Stiles Farm<br>Bulk (n=4)")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL,
       y="Average Relative Abundance (%)",
       title="Class-level 16S rRNA Gene Alpha Diversity") +
  labs(caption = "Other = relative abundance < 1%") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.title = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        plot.caption = element_text(hjust=1.24, size=8),
        plot.title = element_text(hjust=0.5))

ggsave("class.stackedbar.16S.trt.ave.png", width=12, height=8)
